"""Read bounded, staged uploads without extracting untrusted archive paths."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
import zipfile

MAX_STRUCTURES = 100
MAX_FILE_BYTES = 5_000_000
MAX_TOTAL_BYTES = 100_000_000


def _paths(value):
    if value is None or value == "":
        return []
    if isinstance(value, str):
        value = value.strip()
        if value.startswith(("[", '"')):
            try:
                value = json.loads(value)
            except json.JSONDecodeError as exc:
                raise ValueError("The structure upload list is not valid JSON.") from exc
        else:
            return [value]
    if isinstance(value, str):
        return [value]
    if not isinstance(value, list) or not all(isinstance(p, str) for p in value):
        raise ValueError("Structure inputs must be a filename or a list of filenames.")
    return value


def staged_paths(app, environment):
    inputs = (app / "inputs").resolve()
    staged = {}
    for config in (app / "inputs.json", inputs / "inputs.json"):
        if config.is_file():
            staged = json.loads(config.read_text())
            if not isinstance(staged, dict):
                raise ValueError("inputs.json must contain named input values.")
            break
    supplied = []
    for name in ("structures", "structure", "batch_archive"):
        supplied.extend(_paths(environment.get(name) or staged.get(name)))
    if not supplied:
        supplied = [str(p) for p in sorted(inputs.rglob("*"))
                    if p.is_file() and p.name.lower().endswith(
                        (".pdb", ".ent", ".cif", ".mmcif", ".gz", ".tgz", ".zip", ".tar"))]
    resolved = []
    for supplied_path in supplied:
        path = Path(supplied_path)
        if not path.is_absolute():
            path = app / path if path.parts and path.parts[0] == "inputs" else inputs / path
        path = path.resolve()
        if not path.is_relative_to(inputs) or not path.is_file():
            raise ValueError("All uploads must be staged files under /app/inputs/.")
        if path not in resolved:
            resolved.append(path)
    if not resolved:
        raise ValueError("Upload one or more PDB/mmCIF structures or a structure archive.")
    return resolved


def _member_name(name):
    path = PurePosixPath(name.replace("\\", "/"))
    if path.is_absolute() or ".." in path.parts or ":" in str(path):
        raise ValueError("Archive paths must be relative and must not contain traversal.")
    return path.as_posix()


def collect_structures(app, environment, *, max_structures=MAX_STRUCTURES,
                       max_file_bytes=MAX_FILE_BYTES, max_total_bytes=MAX_TOTAL_BYTES):
    structures = []
    total = 0

    def bounded(stream, limit):
        raw = stream.read(limit + 1)
        if len(raw) > limit:
            raise ValueError("An input exceeds the allowed uncompressed size.")
        return raw

    def add(name, raw, provenance=None):
        nonlocal total
        provenance = dict(provenance or {})
        name = _member_name(name)
        if name.lower().endswith(".gz"):
            provenance["compressed_sha256"] = hashlib.sha256(raw).hexdigest()
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
                raw = bounded(stream, max_file_bytes)
            name = name[:-3]
        suffix = Path(name).suffix.lower()
        if suffix not in (".pdb", ".ent", ".cif", ".mmcif"):
            return
        if not raw or len(raw) > max_file_bytes:
            raise ValueError(f"Each structure must be nonempty and at most {max_file_bytes} bytes.")
        if len(structures) >= max_structures:
            raise ValueError(f"A job accepts at most {max_structures} structures.")
        total += len(raw)
        if total > max_total_bytes:
            raise ValueError("The batch exceeds the total uncompressed size limit.")
        try:
            content = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError("Structure files must contain UTF-8 PDB/mmCIF text.") from exc
        if "\x00" in content:
            raise ValueError("Structure files must be text, without NUL bytes.")
        digest = hashlib.sha256(raw).hexdigest()
        if provenance.get("sha256") and provenance["sha256"] != digest:
            raise ValueError("A prepared structure does not match its manifest SHA-256.")
        stem = re.sub(r"[^A-Za-z0-9_-]", "_", Path(name).stem)[:42].strip("_") or "structure"
        key = f"s{len(structures)+1:04d}_{stem}"
        source = str(provenance.get("original_name") or name)
        structures.append({"key": key, "source": source, "structure_text": content,
                           "file_format": "pdb" if suffix in (".pdb", ".ent") else "cif",
                           "input_sha256": digest, "input_bytes": len(raw), "provenance": provenance})

    def read_zip(path):
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            if len(members) > max_structures * 4 + 20:
                raise ValueError("The ZIP contains too many entries.")
            names = set()
            by_name = {}
            for info in members:
                name = _member_name(info.filename)
                if name in names:
                    raise ValueError("The ZIP contains duplicate member names.")
                names.add(name)
                by_name[name] = info
                if stat.S_ISLNK(info.external_attr >> 16):
                    raise ValueError("Symbolic links are not accepted in archives.")
                if info.file_size > max_total_bytes:
                    raise ValueError("A ZIP member exceeds the uncompressed size limit.")
            if "manifest.json" in names:
                with archive.open(by_name["manifest.json"]) as stream:
                    manifest = json.loads(bounded(stream, 1_000_000))
                if not isinstance(manifest, dict):
                    raise ValueError("The prepared archive manifest must be a JSON object.")
                entries = manifest.get("structures")
                if not isinstance(entries, list) or len(entries) > max_structures:
                    raise ValueError("The prepared archive manifest has an invalid structures list.")
                used = set()
                for entry in entries:
                    if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
                        raise ValueError("Each prepared structure needs a relative path.")
                    name = _member_name(entry["path"])
                    if name in used:
                        raise ValueError("A prepared structure occurs twice in the manifest.")
                    used.add(name)
                    if name not in by_name or by_name[name].is_dir():
                        raise ValueError("A prepared structure is missing from the archive.")
                    before = len(structures)
                    with archive.open(by_name[name]) as stream:
                        add(name, bounded(stream, max_file_bytes), entry)
                    if len(structures) == before:
                        raise ValueError("Prepared manifest entries must name supported PDB/mmCIF structures.")
            else:
                for info in members:
                    if info.is_dir():
                        continue
                    name = _member_name(info.filename)
                    if name.lower().endswith((".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".ent.gz", ".cif.gz", ".mmcif.gz")):
                        with archive.open(info) as stream:
                            add(name, bounded(stream, max_file_bytes), {"archive": path.name})

    def read_tar(path):
        with tarfile.open(path, "r:*") as archive:
            count = 0
            seen = set()
            for info in archive:
                count += 1
                if count > max_structures * 4 + 20:
                    raise ValueError("The TAR contains too many entries.")
                name = _member_name(info.name)
                if name in seen:
                    raise ValueError("The TAR contains duplicate member names.")
                seen.add(name)
                if not info.isdir() and not info.isfile():
                    raise ValueError("Only ordinary files and directories are accepted in TAR archives.")
                if info.size > max_total_bytes:
                    raise ValueError("A TAR member exceeds the uncompressed size limit.")
                if info.isfile() and name.lower().endswith((".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".ent.gz", ".cif.gz", ".mmcif.gz")):
                    with archive.extractfile(info) as stream:
                        add(name, bounded(stream, max_file_bytes), {"archive": path.name})

    for path in staged_paths(Path(app), environment):
        if path.stat().st_size > max_total_bytes:
            raise ValueError("The uploaded file exceeds the batch size limit.")
        lower = path.name.lower()
        if lower.endswith(".zip"):
            read_zip(path)
        elif lower.endswith((".tar", ".tar.gz", ".tgz")):
            read_tar(path)
        else:
            with path.open("rb") as stream:
                add(path.name, bounded(stream, max_file_bytes))
    if not structures:
        raise ValueError("No supported PDB/mmCIF structure files were found in the uploads.")
    return structures
