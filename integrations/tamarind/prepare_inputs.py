#!/usr/bin/env python3
"""Prepare Colab-style inputs before uploading them to offline Tamarind jobs.

This companion CLI runs on the user's networked computer, not inside Tamarind.
It accepts public structure IDs, ID-list files, structures, folders, and bounded
archives, then produces a ZIP and a provenance manifest. It never reads an API
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

MAX_STRUCTURES = 100
MAX_FILE_BYTES = 5_000_000
MAX_TOTAL_BYTES = 100_000_000
MAX_ARCHIVE_BYTES = 200_000_000
MAX_ARCHIVE_ENTRIES = 1000
MAX_ID_FILE_BYTES = 1_000_000
PROVIDER_PREFIXES = {
    "files.rcsb.org": "/download/",
    "alphafold.ebi.ac.uk": "/files/",
    "api.esmatlas.com": "/fetchPredictedStructure/",
    "ted.cathdb.info": "/api/v1/files/",
}
UNIPROT = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
AF_MODEL = rf"AF-({UNIPROT})-F([1-9][0-9]*)-model_v([1-9][0-9]*)"
ID_HEADERS = {"id", "ids", "identifier", "identifiers", "structure_id", "structure_ids",
              "query_id", "pdb", "pdb_id", "uniprot", "uniprot_id", "accession"}


def parse_ids(text: str) -> list[str]:
    """Accept comma/whitespace-separated identifiers and optional # comments."""
    clean = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
    return list(dict.fromkeys(part for part in re.split(r"[\s,]+", clean.strip()) if part))


def read_id_list(data: bytes, suffix: str = ".txt") -> list[str]:
    if not data or len(data) > MAX_ID_FILE_BYTES:
        raise ValueError("ID-list files must be nonempty and at most 1 MB.")
    text = data.decode("utf-8-sig")
    if "\x00" in text:
        raise ValueError("An ID list must be a UTF-8 text file.")
    if suffix.lower() == ".txt":
        return parse_ids(text)
    rows = list(csv.reader(io.StringIO(text), delimiter="\t" if suffix.lower() == ".tsv" else ","))
    rows = [row for row in rows if row and any(value.strip() for value in row)
            and not row[0].lstrip().startswith("#")]
    if not rows:
        return []
    header = [value.strip().lower().replace(" ", "_") for value in rows[0]]
    index = next((i for i, value in enumerate(header) if value in ID_HEADERS), None)
    if index is not None:
        values = [row[index] for row in rows[1:] if len(row) > index]
    else:
        values = [value for row in rows for value in row]
    return parse_ids("\n".join(values))


def resolve_identifier(identifier: str) -> dict:
    """Resolve identifiers only through fixed public-provider URL templates."""
    value = identifier.strip()
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


def read_bounded(stream, limit: int, label: str) -> bytes:
    data = stream.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"{label} exceeds its allowed byte limit.")
    return data


def download_structure(url: str, *, max_bytes: int = MAX_FILE_BYTES, attempts: int = 3,
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
                if length is not None and int(length) > max_bytes:
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


def structure_bytes(data: bytes, name: str, limit: int) -> bytes:
    if name.lower().endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as compressed:
            data = read_bounded(compressed, limit, "Decompressed structure")
    if not data or len(data) > limit:
        raise ValueError("Each structure must be nonempty and at most 5 MB.")
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


def archive_structures(path: Path, *, max_file_bytes: int, max_total_bytes: int):
    if path.stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("Input archives must be at most 200 MB.")
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
        if count > MAX_ARCHIVE_ENTRIES or total > max_total_bytes:
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
                    if info.file_size > max_file_bytes:
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
                    if info.size > max_file_bytes:
                        raise ValueError("An archived structure exceeds the allowed byte limit.")
                    with archive.extractfile(info) as stream:
                        yield name, read_bounded(stream, max_file_bytes, "Archived structure")


def is_archive(path: Path) -> bool:
    return path.name.lower().endswith((".zip", ".tar", ".tar.gz", ".tgz"))


def prepare_inputs(output_zip: str | Path, *, ids: str = "", id_files=(), sources=(),
                   downloader=download_structure, max_structures: int = MAX_STRUCTURES,
                   max_file_bytes: int = MAX_FILE_BYTES, max_total_bytes: int = MAX_TOTAL_BYTES,
                   overwrite: bool = False) -> dict:
    """Build a complete prepared ZIP atomically; never write extracted members."""
    if not 1 <= max_structures <= MAX_STRUCTURES or not 1 <= max_file_bytes <= MAX_FILE_BYTES:
        raise ValueError("Preparation limits cannot exceed 100 structures or 5 MB per structure.")
    if not 1 <= max_total_bytes <= MAX_TOTAL_BYTES:
        raise ValueError("Preparation cannot exceed 100 MB of structure data.")
    output = Path(output_zip).resolve()
    manifest_path = output.with_suffix(".manifest.json")
    if output.suffix.lower() != ".zip":
        raise ValueError("The prepared output must use a .zip extension.")
    if not overwrite and (output.exists() or manifest_path.exists()):
        raise FileExistsError("An output already exists; choose another path or use --force.")
    excluded = {output, manifest_path}
    if len(ids.encode("utf-8")) > MAX_ID_FILE_BYTES:
        raise ValueError("Identifier text must be at most 1 MB.")
    identifiers = parse_ids(ids)
    local_entries = []
    collected_bytes = 0

    def collect(name, data, metadata):
        nonlocal collected_bytes
        collected_bytes += len(data)
        if len(local_entries) >= max_structures or collected_bytes > max_total_bytes:
            raise ValueError("Local inputs exceed the configured structure-count or total-byte limit.")
        local_entries.append((name, data, metadata))

    for filename in id_files:
        path = Path(filename)
        if path.is_symlink() or not path.is_file():
            raise ValueError("ID-list inputs must be regular files, not links.")
        with path.open("rb") as stream:
            identifiers.extend(read_id_list(read_bounded(stream, MAX_ID_FILE_BYTES, "ID list"), path.suffix))
    files_seen = 0
    for source in sources:
        path = Path(source)
        if path.is_symlink() or not path.exists():
            raise ValueError("Local inputs must exist and cannot be symbolic links.")
        if path.is_dir():
            paths = []
            for item in path.rglob("*"):
                files_seen += 1
                if files_seen > MAX_ARCHIVE_ENTRIES * 10:
                    raise ValueError("Input folders contain too many entries.")
                paths.append(item)
            paths.sort()
        else:
            paths = [path]
        for item in paths:
            if item.is_symlink():
                raise ValueError("Input folders cannot contain symbolic links.")
            if not item.is_file() or item.resolve() in excluded:
                continue
            if item.suffix.lower() in (".txt", ".csv", ".tsv") and not path.is_dir():
                with item.open("rb") as stream:
                    identifiers.extend(read_id_list(read_bounded(stream, MAX_ID_FILE_BYTES, "ID list"), item.suffix))
            elif is_archive(item):
                for name, data in archive_structures(item, max_file_bytes=max_file_bytes, max_total_bytes=max_total_bytes):
                    collect(name, data, {"source_type": "archive", "archive_name": item.name})
            elif structure_format(item.name):
                with item.open("rb") as stream:
                    collect(item.name, read_bounded(stream, max_file_bytes, "Structure"), {"source_type": "local_file"})
            elif not path.is_dir():
                raise ValueError(f"Unsupported input file type: {item.name}")
    identifiers = list(dict.fromkeys(identifiers))
    resolved = []
    seen_urls = set()
    for identifier in identifiers:
        metadata = resolve_identifier(identifier)
        if metadata["source_url"] not in seen_urls:
            resolved.append((identifier, metadata))
            seen_urls.add(metadata["source_url"])
    if len(resolved) + len(local_entries) > max_structures:
        raise ValueError("Input exceeds the configured structure-count limit.")
    if not resolved and not local_entries:
        raise ValueError("Provide at least one supported identifier or structure.")
    entries, blobs, total_bytes = [], {}, 0

    def add(name: str, data: bytes, metadata: dict):
        nonlocal total_bytes
        data = structure_bytes(data, name, max_file_bytes)
        total_bytes += len(data)
        if total_bytes > max_total_bytes:
            raise ValueError("Combined structure data exceeds the configured total-byte limit.")
        fmt = structure_format(name)
        original = PurePosixPath(name).name
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", original[:-3] if original.lower().endswith(".gz") else original)
        stem = Path(stem).stem.strip("._-")[:100] or "structure"
        archive_path = f"structures/{len(entries) + 1:04d}_{stem}.{fmt}"
        entry = {"path": archive_path, "original_name": original, "file_format": fmt,
                 "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), **metadata}
        if name != original:
            entry["archive_member"] = name
        entries.append(entry)
        blobs[archive_path] = data

    for identifier, metadata in resolved:
        add(metadata["filename"], downloader(metadata["source_url"], max_bytes=max_file_bytes),
            {"source_type": "identifier", "input_id": identifier,
             "canonical_id": metadata["canonical_id"], "source_url": metadata["source_url"]})
    for name, data, metadata in local_entries:
        add(name, data, metadata)
    manifest = {"schema_version": 1, "prepared_by": "ZincSight Tamarind input preparer",
                "requires_network_at_runtime": False, "requested_ids": identifiers,
                "structure_count": len(entries), "total_structure_bytes": total_bytes,
                "structures": entries}
    manifest_bytes = (json.dumps(manifest, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_paths = []
    try:
        for destination in (output, manifest_path):
            with tempfile.NamedTemporaryFile(prefix=".zincsight-", dir=output.parent, delete=False) as stream:
                temporary_paths.append(Path(stream.name))
        with zipfile.ZipFile(temporary_paths[0], "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for name, data in [("manifest.json", manifest_bytes), *blobs.items()]:
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data)
        if temporary_paths[0].stat().st_size > MAX_ARCHIVE_BYTES:
            raise ValueError("Prepared ZIP exceeds the 200 MB upload limit.")
        temporary_paths[1].write_bytes(manifest_bytes)
        os.replace(temporary_paths[0], output)
        os.replace(temporary_paths[1], manifest_path)
    finally:
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="*", help="Structure files, folders, archives, or ID-list files")
    parser.add_argument("--ids", default="", help="Comma/whitespace-separated public structure identifiers")
    parser.add_argument("--ids-file", action="append", default=[], help="UTF-8 .txt/.csv/.tsv ID list (repeatable)")
    parser.add_argument("--output", required=True, type=Path, help="Prepared ZIP to upload as Tamarind's batch archive")
    parser.add_argument("--force", action="store_true", help="Replace an existing prepared output")
    args = parser.parse_args(argv)
    try:
        manifest = prepare_inputs(args.output, ids=args.ids, id_files=args.ids_file,
                                  sources=args.sources, overwrite=args.force)
    except (ValueError, OSError, UnicodeError, zipfile.BadZipFile, tarfile.TarError) as exc:
        parser.exit(2, f"Preparation failed: {exc}\n")
    print(f"Prepared {manifest['structure_count']} structures: {args.output}")
    print(f"Manifest: {args.output.with_suffix('.manifest.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
