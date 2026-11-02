"""Publish a compact Tamarind result view while retaining every detailed output.

The ZIP keeps the original directory layout and JSON paths unchanged. Public
viewer/session names and small tables are staged separately before replacing
individual children, so the output directory remains compatible with mounts.
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
import zipfile

BUNDLE_NAME = "ZincSight_results.zip"
SUMMARY_NAME = "00_Summary.csv"
SITES_NAME = "01_Predicted_sites.csv"
SUMMARY_FIELDS = ("structure", "sites", "top_probability", "status", "structure_file")
SITE_FIELDS = ("structure", "site", "probability", "coordinating_residues",
               "template_id", "template_rmsd_angstrom")
_OLD_BUNDLES = {"zincsight_results.zip", BUNDLE_NAME}


def _member(output, value):
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("Result artifact paths must be relative POSIX filenames.")
    relative = PurePosixPath(value)
    if relative.is_absolute() or ".." in relative.parts or ":" in value:
        raise ValueError("Result artifact paths must remain inside the output directory.")
    path = output / relative
    if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(output):
        raise ValueError(f"A required result artifact is missing or unsafe: {value}")
    return relative.as_posix(), path


def _label(result):
    source = result.get("source")
    if not isinstance(source, str) or not source:
        raise ValueError("Every structure needs its original source name.")
    name = PurePosixPath(source.replace("\\", "/")).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for extension in (".mmcif", ".pdb", ".cif", ".ent"):
        if name.lower().endswith(extension):
            name = name[:-len(extension)]
            break
    # Tamarind prefixes direct upload names with an index. Archive/prepared
    # provenance identifies genuine names, whose leading numbers must remain.
    if not (result.get("input_provenance") or result.get("provenance")):
        name = re.sub(r"^\d+_", "", name, count=1)
    display = name.strip() or "structure"
    filename = re.sub(r"[^A-Za-z0-9._-]+", "_", display).strip("._-")[:64] or "structure"
    return display, filename


def _cell(value):
    if isinstance(value, str) and (value.lstrip().startswith(("=", "+", "-", "@"))
                                   or value.startswith(("\t", "\r", "\n"))):
        return "'" + value
    return value


def _write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: _cell(value) for key, value in row.items()} for row in rows)


def finalize_outputs(output: Path, *, structures: list[dict]) -> dict:
    """Return public filenames and original ZIP members, keyed by structure key.

    ``prediction.json`` is copied unchanged to the public directory. The caller
    can then add public links to that copy; the ZIP retains its original JSON
    and artifact paths for reproducibility.
    """
    output = Path(output)
    if output.is_symlink() or not output.is_dir():
        raise ValueError("The result output must be a regular directory.")
    output = output.resolve()
    if not isinstance(structures, list) or not structures:
        raise ValueError("At least one structure is required to package results.")
    # Capture the original children before creating any temporary files inside
    # the mounted output directory. Never rename the directory itself.
    original_children = sorted(output.iterdir())
    originals = []
    for path in sorted(output.rglob("*")):
        if path.is_symlink():
            raise ValueError("Symbolic links are not accepted in result outputs.")
        if path.is_file() and path.relative_to(output).as_posix() not in _OLD_BUNDLES:
            originals.append((path.relative_to(output).as_posix(), path))
    for required in ("CITATION.txt", "prediction.json", "structures.csv", "sites.csv"):
        _member(output, required)

    mappings, summary_rows, site_rows, public_sources = {}, [], [], []
    width = max(2, len(str(len(structures))))
    for number, result in enumerate(structures, 1):
        if not isinstance(result, dict):
            raise ValueError("Every structure result must be an object.")
        key = result.get("key")
        if (not isinstance(key, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,99}", key)
                or ".." in key or key in mappings):
            raise ValueError("Structure keys must be unique safe identifiers.")
        artifacts = result.get("artifacts")
        if not isinstance(artifacts, dict):
            raise ValueError("Every structure needs its detailed artifact paths.")
        for required in ("mmcif", "pymol_session", "site_annotations"):
            if required not in artifacts:
                raise ValueError(f"Every structure needs a {required} artifact.")
        members = {name: _member(output, path)[0] for name, path in artifacts.items()}
        members["prediction_json"] = _member(output, result.get("prediction_json"))[0]
        members["sites_csv"] = _member(output, f"sites/{key}.csv")[0]
        viewer_kind = "pdb" if "pdb" in artifacts else "mmcif"
        display, label = _label(result)
        base = f"{number:0{width}d}_{label}_ZincSight"
        viewer = base + (".pdb" if viewer_kind == "pdb" else ".cif")
        session = base + ".pse"
        public_sources.extend(((viewer, output / members[viewer_kind]),
                               (session, output / members["pymol_session"])))
        mappings[key] = {"viewer_file": viewer, "pymol_session": session,
                         "bundle_members": members, "display_name": display}
        sites = result.get("sites")
        if not isinstance(sites, list) or not all(isinstance(site, dict) for site in sites):
            raise ValueError("Every structure needs a list of predicted sites.")
        site_count = result.get("site_count")
        if (not isinstance(site_count, int) or isinstance(site_count, bool)
                or site_count < len(sites)):
            raise ValueError("The qualifying site count must include every returned site.")
        for site in sites:
            probability = site.get("probability")
            if (not isinstance(probability, (int, float)) or isinstance(probability, bool)
                    or not math.isfinite(probability) or not 0 <= probability <= 1):
                raise ValueError("Site probabilities must be finite numbers between zero and one.")
            residues = site.get("coordinating_residues")
            if not isinstance(residues, list) or not all(isinstance(residue, str) for residue in residues):
                raise ValueError("Coordinating residues must be a list of identifiers.")
        status = "No sites predicted" if not sites else "Truncated" if result.get("truncated") else "Complete"
        summary_rows.append({"structure": display, "sites": site_count,
                             "top_probability": max((site["probability"] for site in sites), default=""),
                             "status": status, "structure_file": viewer})
        for rank, site in enumerate(sites, 1):
            site_rows.append({"structure": display, "site": rank, "probability": site["probability"],
                              "coordinating_residues": ", ".join(site["coordinating_residues"]),
                              "template_id": site.get("template_id", ""),
                              "template_rmsd_angstrom": site.get("template_rmsd_angstrom", "")})

    metadata = {"bundle_path": BUNDLE_NAME, "structures": mappings,
                "summary_file": SUMMARY_NAME, "predicted_sites_file": SITES_NAME}
    new_names = {SUMMARY_NAME, SITES_NAME, "public_manifest.json"} | {name for name, _ in public_sources}
    if new_names & {name for name, _ in originals}:
        raise ValueError("The output already contains a reserved public result filename.")
    temporary = Path(tempfile.mkdtemp(prefix=".packaging-", dir=output))
    staged, backup = temporary / "public", temporary / "original"
    moved_originals, published = [], []
    rollback_failed = False
    try:
        staged.mkdir()
        backup.mkdir()
        for name in ("CITATION.txt", "prediction.json"):
            shutil.copy2(output / name, staged / name)
        for name, path in public_sources:
            shutil.copy2(path, staged / name)
        _write_csv(staged / SUMMARY_NAME, SUMMARY_FIELDS, summary_rows)
        _write_csv(staged / SITES_NAME, SITE_FIELDS, site_rows)
        with zipfile.ZipFile(staged / BUNDLE_NAME, "w", zipfile.ZIP_DEFLATED) as archive:
            for name, path in originals:
                archive.write(path, name)
            # Public tables/files are convenient after downloading the bundle;
            # the original detailed paths above remain valid and unchanged.
            for name in (SUMMARY_NAME, SITES_NAME, *(name for name, _ in public_sources)):
                archive.write(staged / name, name)
            archive.writestr("public_manifest.json", json.dumps(metadata, indent=2, allow_nan=False) + "\n")
        # /app/out may be a mount or have a parent owned by another user. Move
        # only its children within the same writable volume, retaining every
        # original until all public files have been published.
        try:
            for child in original_children:
                destination = backup / child.name
                os.replace(child, destination)
                moved_originals.append((child, destination))
            for child in sorted(staged.iterdir()):
                destination = output / child.name
                os.replace(child, destination)
                published.append(destination)
        except BaseException as error:
            # Preserve the backup even if recovery itself is interrupted by a
            # non-filesystem exception before it can finish recording errors.
            rollback_failed = True
            rollback_errors = []
            for path in reversed(published):
                try:
                    path.unlink()
                except OSError as exc:
                    rollback_errors.append(exc)
            for original, saved in reversed(moved_originals):
                try:
                    os.replace(saved, original)
                except OSError as exc:
                    rollback_errors.append(exc)
            if rollback_errors:
                raise RuntimeError(
                    f"Result publication failed and rollback was incomplete; "
                    f"original results remain under {output} and {backup}."
                ) from error
            rollback_failed = False
            raise
    finally:
        # A failed rollback must retain the original tree for recovery rather
        # than deleting the only copy while handling the initial exception.
        if not rollback_failed:
            shutil.rmtree(temporary, ignore_errors=True)
    return metadata
