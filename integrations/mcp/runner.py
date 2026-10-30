"""Bounded, isolated calls to the existing ZincSight pipeline (Linux/macOS)."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
MAX_BYTES = 2_000_000
TIMEOUT = 180
LIMITATIONS = (
    "Computational predictions, not experimental confirmation of metal binding. "
    "Scores do not establish metal selectivity, occupancy, or catalytic activity. "
    "No predicted sites does not rule out binding."
)


def get_structure(pdb_id: str | None, structure_text: str | None, file_format: str) -> str:
    if bool(pdb_id) == bool(structure_text):
        raise ValueError("Provide exactly one of pdb_id or structure_text.")
    if file_format not in ("pdb", "cif"):
        raise ValueError("file_format must be pdb or cif.")
    if pdb_id:
        if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
            raise ValueError("pdb_id must be a four-character PDB accession, e.g. 1A1F.")
        if file_format != "pdb":
            raise ValueError("Use file_format=pdb when supplying pdb_id.")
        # Fixed provider; callers cannot supply arbitrary URLs or local paths.
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                content = response.read(MAX_BYTES + 1)
        except Exception as exc:
            raise ValueError("Could not retrieve this PDB structure from RCSB.") from exc
        if len(content) > MAX_BYTES:
            raise ValueError("Structure exceeds the 2 MB limit.")
        structure_text = content.decode("utf-8")
    if len(structure_text.encode("utf-8")) > MAX_BYTES:
        raise ValueError("Structure exceeds the 2 MB limit.")
    if "\x00" in structure_text:
        raise ValueError("Structure must be a text PDB or mmCIF file.")
    if file_format == "pdb" and not any(line.startswith("ATOM  ") for line in structure_text.splitlines()):
        raise ValueError("PDB input must contain protein ATOM records.")
    if file_format == "cif" and "_atom_site." not in structure_text:
        raise ValueError("mmCIF input must contain an atom_site table.")
    return structure_text


def parse_sites(path: Path, minimum_probability: float) -> list[dict]:
    sites = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            probability = float(row["prob"])
            coordinates = json.loads(row["predicted_ion_pos"])
            if not (math.isfinite(probability) and 0 <= probability <= 1):
                raise ValueError("Invalid probability in predictor output.")
            if len(coordinates) != 3 or not all(isinstance(x, (int, float)) and math.isfinite(x) for x in coordinates):
                raise ValueError("Invalid coordinates in predictor output.")
            if probability >= minimum_probability:
                sites.append({
                    "probability": probability,
                    "coordinates_angstrom": coordinates,
                    "coordinating_residues": [row[f"chain_resi_{i}"] for i in range(1, 5) if row[f"chain_resi_{i}"]],
                    "residue_combination": row["resi_comb_ordered"],
                    "template_id": row["template_id"],
                    "template_rmsd_angstrom": float(row["template_rmsd"]),
                })
    return sorted(sites, key=lambda site: site["probability"], reverse=True)


def run_prediction(*, pdb_id: str | None = None, structure_text: str | None = None,
                   file_format: str = "pdb", histidine_rotamers: bool = False,
                   minimum_probability: float = 0.0, max_sites: int = 50) -> dict:
    if not math.isfinite(minimum_probability) or not 0 <= minimum_probability <= 1:
        raise ValueError("minimum_probability must be between 0 and 1.")
    if isinstance(max_sites, bool) or not isinstance(max_sites, int) or not 1 <= max_sites <= 100:
        raise ValueError("max_sites must be an integer between 1 and 100.")
    content = get_structure(pdb_id, structure_text, file_format)
    with tempfile.TemporaryDirectory(prefix="zincsight-mcp-") as directory:
        work = Path(directory)
        # Only scientific source/assets; never copy existing results, DBs, or .env.
        shutil.copytree(ROOT / "src", work / "src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        shutil.copy2(ROOT / "main_execute.py", work / "main_execute.py")
        shutil.copy2(Path(__file__).with_name("worker.py"), work / "worker.py")
        (work / f"input.{file_format}").write_text(content, encoding="utf-8")
        environment = os.environ.copy()
        environment.update(PYTHON_DOTENV_DISABLED="1", DEBUGGING_TEMP="false",
                           KEEP_TEMP_TABLES_TEMP="false", FILTER_PROB_AFTER_COMPRESSION="false",
                           EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD="0.0", OMP_NUM_THREADS="1",
                           OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
        with (work / "worker.log").open("w") as log:
            process = subprocess.Popen(
                [sys.executable, str(work / "worker.py"), file_format, str(int(histidine_rotamers))],
                cwd=work, env=environment, stdout=log, stderr=log, start_new_session=True,
            )
            try:
                returncode = process.wait(timeout=TIMEOUT)
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError("Prediction exceeded the 180-second limit; try a smaller structure.") from exc
            finally:
                # Also stop descendant PostgreSQL/multiprocessing workers before cleanup.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        if returncode != 0 or not (work / "status.json").exists():
            raise RuntimeError("ZincSight failed to process this structure. Verify the input and server dependencies.")
        status = json.loads((work / "status.json").read_text())
        sites = parse_sites(work / "results/table/table_predicted_zn_binding_sites.csv", minimum_probability) if status["has_sites"] else []
        return {
            "method": "ZincSight", "source": pdb_id.upper() if pdb_id else "provided_structure",
            "histidine_rotamers": histidine_rotamers, "minimum_probability": minimum_probability,
            "site_count": len(sites), "sites": sites[:max_sites], "truncated": len(sites) > max_sites,
            "model_policy": "First structural model only; all protein chains in that model.",
            "limitations": LIMITATIONS, "software_url": "https://github.com/MECHTI1/ZincSight",
        }
