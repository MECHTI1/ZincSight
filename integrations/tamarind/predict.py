"""Adapt Tamarind's staged inputs to the pinned ZincSight preview runner."""
from __future__ import annotations

import csv
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace

ENGINE_ROOT = Path("/opt/zincsight")
ENGINE_REVISION = "7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03"
MAX_BYTES = 5_000_000
VALIDATION_STATUS = (
    "Experimental integration with reproducibility fixes; historical outputs can differ. "
    "Published accuracy benchmarks have not been re-established for this revision."
)
CITATION = (
    "Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion "
    "location in proteins with ZincSight. Protein Science. 2025;34(11):e70350. "
    "https://doi.org/10.1002/pro.70350"
)
CSV_FIELDS = [
    "site_rank", "probability", "x_angstrom", "y_angstrom", "z_angstrom",
    "coordinating_residues", "residue_combination", "template_id", "template_rmsd_angstrom",
]


def read_parameters(environment):
    flag = environment.get("histidine_rotamers", "true").strip().lower()
    if flag not in ("true", "false", "1", "0"):
        raise ValueError("histidine_rotamers must be true or false.")
    try:
        probability = float(environment.get("minimum_probability", "0.0001"))
        count = float(environment.get("max_sites", "100"))
    except (TypeError, ValueError) as exc:
        raise ValueError("minimum_probability and max_sites must be numbers.") from exc
    if not math.isfinite(probability) or not 0 <= probability <= 1:
        raise ValueError("minimum_probability must be between 0 and 1.")
    if not math.isfinite(count) or not count.is_integer() or not 1 <= count <= 100:
        raise ValueError("max_sites must be a whole number between 1 and 100.")
    return {"histidine_rotamers": flag in ("true", "1"),
            "minimum_probability": probability, "max_sites": int(count)}


def find_structure(app: Path, environment) -> Path:
    inputs = (app / "inputs").resolve()
    supplied = environment.get("structure", "").strip()
    if supplied:
        path = Path(supplied)
        if not path.is_absolute():
            path = app / path if path.parts[0] == "inputs" else inputs / path
        candidates = [path]
    else:
        # Some Tamarind upload routes only stage the file, without a path env var.
        candidates = [p for p in inputs.rglob("*")
                      if p.is_file() and p.suffix.lower() in (".pdb", ".cif", ".mmcif")]
    if len(candidates) != 1:
        raise ValueError("Upload exactly one PDB or mmCIF structure per job.")
    path = candidates[0].resolve()
    if not path.is_relative_to(inputs) or not path.is_file():
        raise ValueError("The structure must be a staged file under /app/inputs/.")
    if path.suffix.lower() not in (".pdb", ".cif", ".mmcif"):
        raise ValueError("The structure must have a .pdb, .cif or .mmcif extension.")
    return path


@contextmanager
def diagnostic_workdir(*args, **kwargs):
    """Keep job-local worker diagnostics visible before temporary cleanup."""
    with tempfile.TemporaryDirectory(*args, **kwargs) as directory:
        try:
            yield directory
        except Exception:
            for name in ("worker.log", "database.log"):
                path = Path(directory) / name
                if path.is_file():
                    print(f"--- ZincSight {name} (last 32 KiB) ---", file=sys.stderr)
                    print(path.read_text(errors="replace")[-32768:], file=sys.stderr)
            raise


def load_predictor():
    revision_file = ENGINE_ROOT / "SOURCE_REVISION"
    if not revision_file.is_file() or revision_file.read_text().strip() != ENGINE_REVISION:
        raise RuntimeError("The pinned ZincSight engine is missing or mismatched; rebuild the image.")
    # Use the image's pinned scientific source, independently of the mounted repo.
    sys.path.insert(0, str(ENGINE_ROOT))
    from integrations.mcp import runner
    runner.tempfile = SimpleNamespace(TemporaryDirectory=diagnostic_workdir)
    return runner.run_prediction


def run(app: Path, environment, predictor) -> dict:
    parameters = read_parameters(environment)
    structure = find_structure(app, environment)
    with structure.open("rb") as stream:
        raw = stream.read(MAX_BYTES + 1)
    if not raw or len(raw) > MAX_BYTES:
        raise ValueError("The input must be a nonempty structure file of at most 5 MB.")
    try:
        content = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("Upload a text PDB/mmCIF file, without compression.") from exc
    result = predictor(structure_text=content,
                       file_format="pdb" if structure.suffix.lower() == ".pdb" else "cif",
                       database_backend="pglite", **parameters)
    result.update(
        source=structure.name, input_sha256=hashlib.sha256(raw).hexdigest(),
        software_revision=ENGINE_REVISION, integration="tamarind", integration_version="0.1.2",
        database_transport="unix_socket",
        validation_status=VALIDATION_STATUS, citation=CITATION,
        returned_site_count=len(result["sites"]), max_sites=parameters["max_sites"],
    )
    output = app / "out"
    output.mkdir(parents=True, exist_ok=True)
    with (output / "sites.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for rank, site in enumerate(result["sites"], 1):
            x, y, z = site["coordinates_angstrom"]
            writer.writerow({
                "site_rank": rank, "probability": site["probability"],
                "x_angstrom": x, "y_angstrom": y, "z_angstrom": z,
                "coordinating_residues": json.dumps(site["coordinating_residues"]),
                "residue_combination": site["residue_combination"],
                "template_id": site["template_id"],
                "template_rmsd_angstrom": site["template_rmsd_angstrom"],
            })
    (output / "prediction.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Returned {result['returned_site_count']} of {result['site_count']} sites "
          f"passing the probability cutoff; truncated={result['truncated']}.")
    return result


if __name__ == "__main__":
    try:
        run(Path("/app"), os.environ, load_predictor())
    except Exception as exc:
        print(f"ZincSight failed: {exc}", file=sys.stderr)
        sys.exit(1)
