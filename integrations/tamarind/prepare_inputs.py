#!/usr/bin/env python3
"""Prepare Colab-style inputs before uploading them to offline Tamarind jobs.

This companion CLI runs on the user's networked computer, not inside Tamarind.
It accepts public structure IDs, ID-list files, structures, folders, and archives,
then produces upload-sized ZIP parts and a provenance manifest. It never reads an API
key, submits a job, or includes credentials or absolute local paths in outputs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib

UPLOAD_MAX_BYTES = 200_000_000  # Tamarind's per-file upload transport limit.
PROVIDER_PREFIXES = {
    "files.rcsb.org": "/download/",
    "alphafold.ebi.ac.uk": "/files/",
    "api.esmatlas.com": "/fetchPredictedStructure/",
    "ted.cathdb.info": "/api/v1/files/",
}
UNIPROT = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
AF_MODEL = rf"AF-({UNIPROT})-F([1-9][0-9]*)-model_v([1-9][0-9]*)"
ID_HEADERS = {"id", "ids", "identifier", "identifiers", "structure_id", "structure_ids",
              "query_id", "pdb", "pdb_id", "uniprot", "uniprot_id", "accession",
              "entry", "uniprot_accession", "model_id", "alphafold_model_id"}


def _id_text(text: str) -> str:
    if not isinstance(text, str) or re.search(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", text):
        raise ValueError("ID lists must contain text without binary control characters.")
    return text.strip().lstrip("\ufeff").strip()


def _unquote_id(value: str) -> str:
    value = _id_text(value)
    if value.startswith(('"', "'")):
        if len(value) < 2 or value[-1] != value[0]:
            raise ValueError("An identifier has unmatched quotes.")
        value = value[1:-1].strip().lstrip("\ufeff").strip()
    if '"' in value or "'" in value:
        raise ValueError("An identifier has malformed quoting.")
    return value


def parse_ids(text: str) -> list[str]:
    """Accept quoted/BOM-prefixed comma/whitespace IDs and optional # comments."""
    clean = "\n".join(line.split("#", 1)[0] for line in _id_text(text).splitlines())
    tokens, previous_end = [], None
    for match in re.finditer(r'''(?:"[^"]*"|'[^']*'|[^\s,]+)''', clean):
        if previous_end is not None and not clean[previous_end:match.start()]:
            raise ValueError("Quoted identifiers must be separated by commas or whitespace.")
        value = _unquote_id(match.group())
        tokens.extend(part for part in re.split(r"[\s,]+", value) if part)
        previous_end = match.end()
    return list(dict.fromkeys(tokens))


def read_id_list(data: bytes, suffix: str = ".txt", *, max_bytes=None) -> list[str]:
    if not data or (max_bytes is not None and len(data) > max_bytes):
        raise ValueError("ID-list files must be nonempty and respect any configured byte guard.")
    text = _id_text(data.decode("utf-8-sig"))
    # Detect the table from its header, not only its filename: spreadsheets are
    # often saved as .txt or given the wrong .csv/.tsv extension. Keep the first
    # recognized ID column; metadata columns never become extra identifiers.
    delimiters = ("\t", ",") if suffix.lower() == ".tsv" else (",", "\t")
    for delimiter in delimiters:
        rows = (row for row in csv.reader(io.StringIO(text), delimiter=delimiter,
                                           skipinitialspace=True, strict=True)
                if row and any(value.strip() for value in row)
                and not row[0].lstrip().startswith("#"))
        try:
            first = next(rows, None)
        except csv.Error:
            continue
        if first is None:
            continue
        header = [re.sub(r"\s+", "_", value.strip().lower()) for value in first]
        index = next((i for i, value in enumerate(header) if value in ID_HEADERS), None)
        if index is not None:
            try:
                values = []
                for row in rows:
                    selected = parse_ids(row[index]) if len(row) > index else []
                    if not selected:
                        raise ValueError(f"A table row is missing a value for the selected ID column: {header[index]}.")
                    values.extend(selected)
            except csv.Error as exc:
                raise ValueError("The ID-list table has malformed quoting.") from exc
            return list(dict.fromkeys(values))
    return parse_ids(text)


def resolve_identifier(identifier: str) -> dict:
    """Resolve identifiers only through fixed public-provider URL templates."""
    value = _unquote_id(identifier)
    value = re.sub(r"\.(?:pdb|cif|mmcif|ent)(?:\.gz)?$", "", value, flags=re.I)
    pdb = re.fullmatch(r"([0-9][A-Z0-9]{3})(?:-assembly([1-9][0-9]{0,2}))?", value, re.I)
    if pdb:
        canonical = pdb[1].upper() + (f"-assembly{pdb[2]}" if pdb[2] else "")
        return {"canonical_id": canonical, "source_url": f"https://files.rcsb.org/download/{canonical}.cif",
                "file_format": "cif", "filename": canonical + ".cif"}
    if re.fullmatch(r"MGYP[0-9]{12}", value, re.I):
        canonical = value.upper()
        return {"canonical_id": canonical,
                "source_url": f"https://api.esmatlas.com/fetchPredictedStructure/{canonical}.pdb",
                "file_format": "pdb", "filename": canonical + ".pdb"}
    ted = re.fullmatch(AF_MODEL + r"_TED([0-9]{2,4})", value, re.I)
    if ted and int(ted[4]) > 0:
        canonical = f"AF-{ted[1].upper()}-F{ted[2]}-model_v{ted[3]}_TED{ted[4]}"
        return {"canonical_id": canonical, "source_url": f"https://ted.cathdb.info/api/v1/files/{canonical}.pdb",
                "file_format": "pdb", "filename": canonical + ".pdb"}
    model = re.fullmatch(AF_MODEL, value, re.I)
    if model:
        canonical = f"AF-{model[1].upper()}-F{model[2]}-model_v{model[3]}"
    elif re.fullmatch(UNIPROT, value, re.I):
        canonical = f"AF-{value.upper()}-F1-model_v6"
    else:
        raise ValueError("Unsupported identifier; use a PDB/assembly, UniProt, full AlphaFold, MGYP or TED domain ID.")
    return {"canonical_id": canonical, "source_url": f"https://alphafold.ebi.ac.uk/files/{canonical}.cif",
            "file_format": "cif", "filename": canonical + ".cif"}


def validate_provider_url(url: str) -> None:
    parsed = urllib.parse.urlsplit(url)
    prefix = PROVIDER_PREFIXES.get(parsed.hostname or "")
    if (parsed.scheme != "https" or prefix is None or not parsed.path.startswith(prefix)
            or parsed.username is not None or parsed.password is not None
            or parsed.port not in (None, 443) or parsed.query or parsed.fragment):
        raise ValueError("Only fixed HTTPS structure-provider URLs are permitted.")


class ProviderRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, msg, headers, newurl):
        validate_provider_url(newurl)
        return super().redirect_request(request, fp, code, msg, headers, newurl)


def read_bounded(stream, limit: int | None, label: str) -> bytes:
    data = stream.read() if limit is None else stream.read(limit + 1)
    if limit is not None and len(data) > limit:
        raise ValueError(f"{label} exceeds its allowed byte limit.")
    return data


def download_structure(url: str, *, max_bytes: int | None = None, attempts: int = 3,
                       timeout: int = 30, opener=None, sleep=time.sleep) -> bytes:
    validate_provider_url(url)
    if not 1 <= attempts <= 3 or not 1 <= timeout <= 60:
        raise ValueError("Downloads allow 1–3 attempts and a 1–60 second timeout.")
    opener = opener or urllib.request.build_opener(ProviderRedirectHandler())
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "ZincSight-Tamarind-Input-Preparer/1.0"})
            with opener.open(request, timeout=timeout) as response:
                validate_provider_url(response.geturl())
                length = response.headers.get("Content-Length")
                if length is not None and max_bytes is not None and int(length) > max_bytes:
                    raise ValueError("Downloaded structure exceeds the allowed byte limit.")
                if "text/html" in response.headers.get("Content-Type", "").lower():
                    raise ValueError("The provider returned an HTML page instead of a structure.")
                return read_bounded(response, max_bytes, "Downloaded structure")
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt + 1 == attempts:
                raise ValueError(f"Public structure provider returned HTTP {exc.code}.") from exc
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt + 1 == attempts:
                raise ValueError("Public structure download failed after bounded retries.") from exc
        sleep(0.25 * (2 ** attempt))
    raise AssertionError("Unreachable download retry state")


def structure_format(name: str) -> str | None:
    lower = name.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    suffix = PurePosixPath(lower).suffix
    return "pdb" if suffix in (".pdb", ".ent") else "cif" if suffix in (".cif", ".mmcif") else None


def structure_bytes(data: bytes, name: str, limit: int | None) -> bytes:
    if name.lower().endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as compressed:
            data = read_bounded(compressed, limit, "Decompressed structure")
    if not data or (limit is not None and len(data) > limit):
        raise ValueError("Each structure must be nonempty and respect any configured byte guard.")
    text = data.decode("utf-8-sig")
    if "\x00" in text:
        raise ValueError("Structures must be text PDB or mmCIF files.")
    fmt = structure_format(name)
    if fmt == "pdb" and not any(line.startswith("ATOM  ") for line in text.splitlines()):
        raise ValueError(f"{PurePosixPath(name).name}: PDB input must contain protein ATOM records.")
    if fmt == "cif" and "_atom_site." not in text:
        raise ValueError(f"{PurePosixPath(name).name}: mmCIF input must contain an atom_site table.")
    return data


def safe_archive_name(name: str) -> str:
    clean = name.rstrip("/")
    # Standard tar tools commonly emit ./member paths; normalize that prefix
    # before detecting traversal and collisions, without allowing internal dots.
    while clean.startswith("./"):
        clean = clean[2:]
    if (not clean or "\\" in clean or clean.startswith("/") or re.match(r"^[A-Za-z]:", clean)
            or any(part in ("", ".", "..") for part in clean.split("/"))):
        raise ValueError("Archive paths must be relative and cannot contain traversal components.")
    return PurePosixPath(clean).as_posix()


def archive_structures(path: Path, *, max_file_bytes=None, max_total_bytes=None, max_entries=None):
    seen, count, total = set(), 0, 0

    def admit(name: str, size: int):
        nonlocal count, total
        name = safe_archive_name(name)
        key = name.casefold()
        if key in seen:
            raise ValueError("Archive contains duplicate or case-colliding paths.")
        seen.add(key)
        count += 1
        total += size
        if ((max_entries is not None and count > max_entries)
                or (max_total_bytes is not None and total > max_total_bytes)):
            raise ValueError("Archive exceeds its entry-count or uncompressed-byte limit.")
        return name

    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if info.is_dir() and info.filename in (".", "./"):
                    continue
                name = admit(info.filename, info.file_size)
                kind = stat.S_IFMT(info.external_attr >> 16)
                if kind not in (0, stat.S_IFREG, stat.S_IFDIR) or info.flag_bits & 1:
                    raise ValueError("Archives cannot contain links, special files or encrypted entries.")
                if info.is_dir():
                    continue
                if structure_format(name):
                    if max_file_bytes is not None and info.file_size > max_file_bytes:
                        raise ValueError("An archived structure exceeds the allowed byte limit.")
                    with archive.open(info) as stream:
                        yield name, read_bounded(stream, max_file_bytes, "Archived structure")
    else:
        with tarfile.open(path, "r|*") as archive:
            for info in archive:
                if info.isdir() and info.name in (".", "./"):
                    continue
                name = admit(info.name, info.size)
                if not (info.isfile() or info.isdir()):
                    raise ValueError("Archives cannot contain links or special files.")
                if info.isfile() and structure_format(name):
                    if max_file_bytes is not None and info.size > max_file_bytes:
                        raise ValueError("An archived structure exceeds the allowed byte limit.")
                    with archive.extractfile(info) as stream:
                        yield name, read_bounded(stream, max_file_bytes, "Archived structure")


def is_archive(path: Path) -> bool:
    return path.name.lower().endswith((".zip", ".tar", ".tar.gz", ".tgz"))


def _guard(value, label):
    if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
        raise ValueError(f"{label} must be a positive integer when supplied.")


def _manifest(entries, identifiers):
    return {"schema_version": 1, "prepared_by": "ZincSight Tamarind input preparer",
            "requires_network_at_runtime": False, "requested_ids": identifiers,
            "structure_count": len(entries),
            "total_structure_bytes": sum(entry["size_bytes"] for entry in entries),
            "structures": entries}


def _json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()


def _compressed_size(path):
    compressor, size = zlib.compressobj(level=6, wbits=-15), 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(compressor.compress(block))
    return size + len(compressor.flush())


def _write_archive(path, manifest, blobs):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, blob in [("manifest.json", _json_bytes(manifest)),
                           *((entry["path"], blobs[entry["path"]]) for entry in manifest["structures"])]:
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type, info.external_attr = zipfile.ZIP_DEFLATED, 0o100644 << 16
            with archive.open(info, "w") as target:
                if isinstance(blob, bytes):
                    target.write(blob)
                else:
                    with blob.open("rb") as stream:
                        for block in iter(lambda: stream.read(1024 * 1024), b""):
                            target.write(block)


def prepare_inputs(output_zip: str | Path, *, ids: str = "", id_files=(), sources=(),
                   downloader=download_structure, max_structures=None, max_file_bytes=None,
                   max_total_bytes=None, max_archive_entries=None, max_folder_entries=None,
                   max_id_file_bytes=None, upload_max_bytes=UPLOAD_MAX_BYTES,
                   overwrite: bool = False) -> dict:
    """Prepare every input, splitting only by the upload transport's byte limit.

    Count, file, total, and archive/folder guards are optional user choices.
    Structure content is spooled to disk rather than retained for the whole batch.
    Every returned archive has its own valid, hash-bearing subset manifest.
    """
    for label, value in (("Structure-count guard", max_structures), ("Structure-byte guard", max_file_bytes),
                         ("Total-byte guard", max_total_bytes), ("Archive-entry guard", max_archive_entries),
                         ("Folder-entry guard", max_folder_entries), ("ID-list byte guard", max_id_file_bytes),
                         ("Upload-byte limit", upload_max_bytes)):
        _guard(value, label)
    if upload_max_bytes is None or upload_max_bytes > UPLOAD_MAX_BYTES:
        raise ValueError("Tamarind upload parts cannot exceed 200 MB.")
    output = Path(output_zip).resolve()
    manifest_path = output.with_suffix(".manifest.json")
    if output.suffix.lower() != ".zip":
        raise ValueError("The prepared output must use a .zip extension.")
    if not overwrite and (output.exists() or manifest_path.exists()):
        raise FileExistsError("An output already exists; choose another path or use --force.")
    if max_id_file_bytes is not None and len(ids.encode()) > max_id_file_bytes:
        raise ValueError("Identifier text exceeds the configured byte guard.")
    identifiers = parse_ids(ids)
    entries, blobs, local_entries, total_bytes = [], {}, [], 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="zincsight-preparation-") as scratch:
        work = Path(scratch)

        def read_ids(path):
            if path.is_symlink() or not path.is_file():
                raise ValueError("ID-list inputs must be regular files, not links.")
            with path.open("rb") as stream:
                return read_id_list(read_bounded(stream, max_id_file_bytes, "ID list"), path.suffix,
                                    max_bytes=max_id_file_bytes)

        def collect(name, data, metadata):
            if max_structures is not None and len(local_entries) >= max_structures:
                raise ValueError("Local inputs exceed the configured structure-count guard.")
            raw_path = work / f"local-{len(local_entries)}"
            raw_path.write_bytes(data)
            local_entries.append((name, raw_path, metadata))

        def is_own_output(path):
            path = path.resolve()
            return (path in (output, manifest_path) or path.is_relative_to(work)
                    or (path.parent == output.parent
                        and re.fullmatch(re.escape(output.stem) + r"\.part[0-9]+\.zip", path.name)))

        for filename in id_files:
            identifiers.extend(read_ids(Path(filename)))
        files_seen = 0
        for source in sources:
            path = Path(source)
            if path.is_symlink() or not path.exists():
                raise ValueError("Local inputs must exist and cannot be symbolic links.")
            if path.is_dir():
                paths = []
                for item in path.rglob("*"):
                    if is_own_output(item):
                        continue
                    files_seen += 1
                    if max_folder_entries is not None and files_seen > max_folder_entries:
                        raise ValueError("Input folders exceed the configured entry-count guard.")
                    paths.append(item)
                paths.sort()
            else:
                paths = [path]
            for item in paths:
                if item.is_symlink():
                    raise ValueError("Input folders cannot contain symbolic links.")
                if not item.is_file() or is_own_output(item):
                    continue
                if item.suffix.lower() in (".txt", ".csv", ".tsv") and not path.is_dir():
                    identifiers.extend(read_ids(item))
                elif is_archive(item):
                    for name, data in archive_structures(item, max_file_bytes=max_file_bytes,
                                                         max_total_bytes=max_total_bytes, max_entries=max_archive_entries):
                        collect(name, data, {"source_type": "archive", "archive_name": item.name})
                elif structure_format(item.name):
                    with item.open("rb") as stream:
                        collect(item.name, read_bounded(stream, max_file_bytes, "Structure"), {"source_type": "local_file"})
                elif not path.is_dir():
                    raise ValueError(f"Unsupported input file type: {item.name}")
        identifiers = list(dict.fromkeys(identifiers))
        resolved, seen_urls = [], set()
        for identifier in identifiers:
            metadata = resolve_identifier(identifier)
            if metadata["source_url"] not in seen_urls:
                resolved.append((identifier, metadata))
                seen_urls.add(metadata["source_url"])
        if max_structures is not None and len(resolved) + len(local_entries) > max_structures:
            raise ValueError("Input exceeds the configured structure-count guard.")
        if not resolved and not local_entries:
            raise ValueError("Provide at least one supported identifier or structure.")

        def add(name, data, metadata):
            nonlocal total_bytes
            data = structure_bytes(data, name, max_file_bytes)
            total_bytes += len(data)
            if max_total_bytes is not None and total_bytes > max_total_bytes:
                raise ValueError("Combined structure data exceeds the configured total-byte guard.")
            fmt, original = structure_format(name), PurePosixPath(name).name
            stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", original[:-3] if original.lower().endswith(".gz") else original)
            stem = Path(stem).stem.strip("._-")[:100] or "structure"
            archive_path = f"structures/{len(entries) + 1:04d}_{stem}.{fmt}"
            entry = {"path": archive_path, "original_name": original, "file_format": fmt,
                     "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), **metadata}
            if name != original:
                entry["archive_member"] = name
            blob = work / f"structure-{len(entries)}"
            blob.write_bytes(data)
            entries.append(entry)
            blobs[archive_path] = blob

        for identifier, metadata in resolved:
            add(metadata["filename"], downloader(metadata["source_url"], max_bytes=max_file_bytes),
                {"source_type": "identifier", "input_id": identifier,
                 "canonical_id": metadata["canonical_id"], "source_url": metadata["source_url"]})
        for name, raw_path, metadata in local_entries:
            add(name, raw_path.read_bytes(), metadata)
            raw_path.unlink()

        # Use measured compressed structure sizes plus a conservative metadata
        # allowance to form upload parts, then verify the actual final ZIP size.
        groups, current, estimate = [], [], 1024
        for entry in entries:
            member_cost = (_compressed_size(blobs[entry["path"]]) + 76 + 2 * len(entry["path"].encode())
                           + len(_json_bytes(entry)) + 256)
            if current and estimate + member_cost > upload_max_bytes:
                groups.append(current)
                current, estimate = [], 1024
            current.append(entry)
            estimate += member_cost
        if current:
            groups.append(current)

        while True:
            archives, manifests, split = [], [], None
            for index, group in enumerate(groups, 1):
                name = output.name if len(groups) == 1 else f"{output.stem}.part{index:03d}.zip"
                for entry in group:
                    entry["archive_file"] = name
                archive_record = {"path": name, "structure_count": len(group),
                                  "total_structure_bytes": sum(entry["size_bytes"] for entry in group)}
                part_ids = list(dict.fromkeys(entry["input_id"] for entry in group if "input_id" in entry))
                manifest = _manifest(group, identifiers if len(groups) == 1 else part_ids)
                manifest["archives"] = [archive_record]
                candidate = work / f"part-{index}.zip"
                _write_archive(candidate, manifest, blobs)
                if candidate.stat().st_size > upload_max_bytes:
                    if len(group) == 1:
                        raise ValueError("A single compressed structure cannot fit the configured Tamarind upload-byte limit; no inputs were omitted.")
                    split = index - 1
                    break
                archives.append((output.parent / name, candidate))
                manifests.append(archive_record)
            if split is None:
                break
            group = groups[split]
            groups[split:split + 1] = [group[:len(group) // 2], group[len(group) // 2:]]

        manifest = _manifest(entries, identifiers)
        manifest["archives"] = manifests
        staged_manifest = work / "manifest.json"
        staged_manifest.write_bytes(_json_bytes(manifest))
        publications = [*archives, (manifest_path, staged_manifest)]
        if not overwrite and any(destination.exists() for destination, _ in publications):
            raise FileExistsError("An output already exists; choose another path or use --force.")
        # Stage on the output filesystem, retaining previous files for rollback.
        with tempfile.TemporaryDirectory(prefix=".zincsight-", dir=output.parent) as transaction:
            staging, backups, published = Path(transaction), [], []
            for index, (_, original) in enumerate(publications):
                with original.open("rb") as source_stream, (staging / str(index)).open("wb") as target:
                    for block in iter(lambda: source_stream.read(1024 * 1024), b""):
                        target.write(block)
            try:
                destinations = {destination for destination, _ in publications}
                obsolete = [candidate for candidate in output.parent.iterdir()
                            if candidate.is_file() and candidate not in destinations
                            and (candidate == output or re.fullmatch(
                                re.escape(output.stem) + r"\.part[0-9]+\.zip", candidate.name))]
                if overwrite:
                    for index, destination in enumerate(obsolete):
                        backup = staging / f"obsolete-{index}"
                        os.replace(destination, backup)
                        backups.append((destination, backup))
                for index, (destination, _) in enumerate(publications):
                    if destination.exists():
                        backup = staging / f"backup-{index}"
                        os.replace(destination, backup)
                        backups.append((destination, backup))
                    os.replace(staging / str(index), destination)
                    published.append(destination)
            except BaseException:
                for destination in reversed(published):
                    destination.unlink(missing_ok=True)
                for destination, backup in reversed(backups):
                    os.replace(backup, destination)
                raise
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="*", help="Structure files, folders, archives, or ID-list files")
    parser.add_argument("--ids", default="", help="Comma/whitespace-separated public structure identifiers")
    parser.add_argument("--ids-file", action="append", default=[], help="UTF-8 .txt/.csv/.tsv ID list; may be repeated")
    parser.add_argument("--output", required=True, type=Path, help="ZIP basename; large batches produce numbered upload parts")
    parser.add_argument("--force", action="store_true", help="Replace an existing prepared output family")
    for option, help_text in (("max-structures", "structure count"), ("max-file-bytes", "bytes per structure"),
                              ("max-total-bytes", "total uncompressed structure bytes"), ("max-archive-entries", "archive entries"),
                              ("max-folder-entries", "folder entries"), ("max-id-file-bytes", "ID-list bytes")):
        parser.add_argument("--" + option, type=int, default=None, help="Optional guard: " + help_text + "; disabled by default")
    parser.add_argument("--upload-max-bytes", type=int, default=UPLOAD_MAX_BYTES,
                        help="Maximum bytes per prepared ZIP (default: Tamarind's 200 MB upload limit)")
    args = parser.parse_args(argv)
    try:
        manifest = prepare_inputs(args.output, ids=args.ids, id_files=args.ids_file, sources=args.sources,
                                  overwrite=args.force, max_structures=args.max_structures,
                                  max_file_bytes=args.max_file_bytes, max_total_bytes=args.max_total_bytes,
                                  max_archive_entries=args.max_archive_entries, max_folder_entries=args.max_folder_entries,
                                  max_id_file_bytes=args.max_id_file_bytes, upload_max_bytes=args.upload_max_bytes)
    except (ValueError, OSError, UnicodeError, zipfile.BadZipFile, tarfile.TarError) as exc:
        parser.exit(2, f"Preparation failed: {exc}\n")
    print(f"Prepared all {manifest['structure_count']} structures in {len(manifest['archives'])} upload archive(s):")
    for archive in manifest["archives"]:
        print(args.output.parent / archive["path"])
    print(f"Manifest: {args.output.with_suffix('.manifest.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
