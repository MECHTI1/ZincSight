"""Stream staged structures and safe archives without scientific batch-size caps."""
from __future__ import annotations

import codecs
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
import zipfile


_FORMATS = ('.pdb', '.ent', '.cif', '.mmcif')


def _paths(value):
    if value is None or value == '':
        return []
    if isinstance(value, str):
        value = value.strip()
        if value.startswith(('[', '"')):
            try:
                value = json.loads(value)
            except json.JSONDecodeError as exc:
                raise ValueError('The structure upload list is not valid JSON.') from exc
        else:
            return [value]
    if isinstance(value, str):
        return [value]
    if not isinstance(value, list) or not all(isinstance(p, str) for p in value):
        raise ValueError('Structure inputs must be a filename or a list of filenames.')
    return value


def staged_paths(app, environment):
    inputs = (app / 'inputs').resolve()
    staged = {}
    for config in (app / 'inputs.json', inputs / 'inputs.json'):
        if config.is_file():
            staged = json.loads(config.read_text())
            if not isinstance(staged, dict):
                raise ValueError('inputs.json must contain named input values.')
            break
    supplied = []
    for name in ('structures', 'structure', 'batch_archive'):
        supplied.extend(_paths(environment.get(name) or staged.get(name)))
    if not supplied:
        supplied = [str(p) for p in sorted(inputs.rglob('*'))
                    if p.is_file() and p.name.lower().endswith(
                        ('.pdb', '.ent', '.cif', '.mmcif', '.gz', '.tgz', '.zip', '.tar'))]
    resolved = []
    for supplied_path in supplied:
        path = Path(supplied_path)
        if not path.is_absolute():
            path = app / path if path.parts and path.parts[0] == 'inputs' else inputs / path
        path = path.resolve()
        if not path.is_relative_to(inputs) or not path.is_file():
            raise ValueError('All uploads must be staged files under /app/inputs/.')
        if path not in resolved:
            resolved.append(path)
    if not resolved:
        raise ValueError('Upload PDB/mmCIF structures or prepared structure archives.')
    return resolved


def _member_name(name):
    path = PurePosixPath(name.replace('\\', '/'))
    if path.is_absolute() or '..' in path.parts or ':' in str(path):
        raise ValueError('Archive paths must be relative and must not contain traversal.')
    return path.as_posix()


def _supported(name):
    name = name.lower()
    return name.endswith(_FORMATS) or name.endswith(tuple(suffix + '.gz' for suffix in _FORMATS))


class _DigestReader:
    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        raw = self.stream.read(size)
        self.digest.update(raw)
        return raw


def collect_structures(app, environment, *, max_structures=None, max_file_bytes=None,
                       max_total_bytes=None, staging_dir=None):
    """Collect every input; production stages bytes on disk, not in a batch string.

    Optional guards are available to library callers. None imposes no artificial
    structure-count, file-size or aggregate-size ceiling. Actual job disk, RAM,
    upload and execution-time resources still apply.
    """
    for name, value in (('max_structures', max_structures), ('max_file_bytes', max_file_bytes),
                        ('max_total_bytes', max_total_bytes)):
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
            raise ValueError(f'{name} must be a positive integer or None.')
    stage = Path(staging_dir) if staging_dir is not None else None
    if stage is not None:
        if stage.is_symlink():
            raise ValueError('Use a private ordinary staging directory.')
        stage.mkdir(parents=True, exist_ok=True)
        stage = stage.resolve()
    structures, total = [], 0

    def add(name, stream, provenance=None):
        nonlocal total
        name = _member_name(name)
        if not _supported(name):
            return
        if max_structures is not None and len(structures) >= max_structures:
            raise ValueError(f'The requested structure-count guard is {max_structures}.')
        provenance = dict(provenance or {})
        compressed = name.lower().endswith('.gz')
        compressed_reader = _DigestReader(stream) if compressed else None
        content = gzip.GzipFile(fileobj=compressed_reader) if compressed else stream
        if compressed:
            name = name[:-3]
        suffix = Path(name).suffix.lower()
        file_format = 'pdb' if suffix in ('.pdb', '.ent') else 'cif'
        stem = re.sub(r'[^A-Za-z0-9_-]', '_', Path(name).stem)[:42].strip('_') or 'structure'
        key = f's{len(structures)+1:04d}_{stem}'
        destination = stage / f'{key}.{file_format}' if stage is not None else None
        sink = destination.open('xb') if destination is not None else io.BytesIO()
        digest, decoder, size = hashlib.sha256(), codecs.getincrementaldecoder('utf-8')(), 0
        try:
            while True:
                raw = content.read(1024 * 1024)
                if not raw:
                    break
                size += len(raw)
                if max_file_bytes is not None and size > max_file_bytes:
                    raise ValueError('An input exceeds the requested uncompressed size guard.')
                if max_total_bytes is not None and total + size > max_total_bytes:
                    raise ValueError('The batch exceeds the requested total uncompressed size guard.')
                if b'\x00' in raw:
                    raise ValueError('Structure files must be text, without NUL bytes.')
                decoder.decode(raw, final=False)
                digest.update(raw)
                sink.write(raw)
            decoder.decode(b'', final=True)
            if not size:
                raise ValueError('Each structure must contain nonempty PDB/mmCIF text.')
            checksum = digest.hexdigest()
            if provenance.get('sha256') and provenance['sha256'] != checksum:
                raise ValueError('A prepared structure does not match its manifest SHA-256.')
            if compressed:
                provenance['compressed_sha256'] = compressed_reader.digest.hexdigest()
            item = {'key': key, 'source': str(provenance.get('original_name') or name),
                    'file_format': file_format, 'input_sha256': checksum,
                    'input_bytes': size, 'provenance': provenance}
            if destination is None:
                item['structure_text'] = sink.getvalue().decode('utf-8')
            else:
                item['structure_path'] = str(destination)
            structures.append(item)
            total += size
        except UnicodeDecodeError as exc:
            if destination is not None:
                destination.unlink(missing_ok=True)
            raise ValueError('Structure files must contain UTF-8 PDB/mmCIF text.') from exc
        except OSError as exc:
            if destination is not None:
                destination.unlink(missing_ok=True)
            if exc.errno == 28:
                raise ValueError('Input staging exhausted available job disk; split this batch across jobs.') from exc
            raise
        except BaseException:
            if destination is not None:
                destination.unlink(missing_ok=True)
            raise
        finally:
            sink.close()
            if compressed:
                content.close()

    def read_zip(path):
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            if max_structures is not None and len(members) > max_structures * 4 + 20:
                raise ValueError('The ZIP contains too many entries for the requested guard.')
            by_name = {}
            for info in members:
                name = _member_name(info.filename)
                if name in by_name:
                    raise ValueError('The ZIP contains duplicate member names.')
                by_name[name] = info
                if stat.S_ISLNK(info.external_attr >> 16):
                    raise ValueError('Symbolic links are not accepted in archives.')
                if info.flag_bits & 1:
                    raise ValueError('Encrypted ZIP members are not supported.')
            if 'manifest.json' in by_name:
                with archive.open(by_name['manifest.json']) as stream:
                    manifest = json.load(io.TextIOWrapper(stream, encoding='utf-8'))
                if not isinstance(manifest, dict):
                    raise ValueError('The prepared archive manifest must be a JSON object.')
                entries = manifest.get('structures')
                if (not isinstance(entries, list) or (max_structures is not None
                        and len(entries) > max_structures)):
                    raise ValueError('The prepared archive manifest has an invalid structures list.')
                used = set()
                for entry in entries:
                    if not isinstance(entry, dict) or not isinstance(entry.get('path'), str):
                        raise ValueError('Each prepared structure needs a relative path.')
                    name = _member_name(entry['path'])
                    if name in used:
                        raise ValueError('A prepared structure occurs twice in the manifest.')
                    used.add(name)
                    if name not in by_name or by_name[name].is_dir():
                        raise ValueError('A prepared structure is missing from the archive.')
                    if not _supported(name):
                        raise ValueError('Prepared manifest entries must name supported PDB/mmCIF structures.')
                    with archive.open(by_name[name]) as stream:
                        add(name, stream, entry)
            else:
                for info in members:
                    if not info.is_dir() and _supported(info.filename):
                        with archive.open(info) as stream:
                            add(info.filename, stream, {'archive': path.name})

    def read_tar(path):
        with tarfile.open(path, 'r:*') as archive:
            seen = set()
            for count, info in enumerate(archive, 1):
                if max_structures is not None and count > max_structures * 4 + 20:
                    raise ValueError('The TAR contains too many entries for the requested guard.')
                name = _member_name(info.name)
                if name in seen:
                    raise ValueError('The TAR contains duplicate member names.')
                seen.add(name)
                if not info.isdir() and not info.isfile():
                    raise ValueError('Only ordinary files and directories are accepted in TAR archives.')
                if info.isfile() and _supported(name):
                    with archive.extractfile(info) as stream:
                        add(name, stream, {'archive': path.name})

    for path in staged_paths(Path(app), environment):
        lower = path.name.lower()
        if lower.endswith('.zip'):
            read_zip(path)
        elif lower.endswith(('.tar', '.tar.gz', '.tgz')):
            read_tar(path)
        else:
            with path.open('rb') as stream:
                add(path.name, stream)
    if not structures:
        raise ValueError('No supported PDB/mmCIF structure files were found in the uploads.')
    return structures
