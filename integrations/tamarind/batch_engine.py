"""Run a bounded structure batch through one isolated, pinned ZincSight engine.

The scientific pipeline already accepts a list of structures. A single private
database per batch avoids paying its startup and template-search costs per file.
"""
from __future__ import annotations

import csv
import importlib.util
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

MAX_STRUCTURES = 100
MAX_BYTES = 5_000_000
MAX_TOTAL_BYTES = 100_000_000
MAX_RESIDUES = 4000
MODEL_POLICY = "First structural model only; all protein chains in that model."
KEY_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z")
RESULT_TABLE = Path("results/table/table_predicted_zn_binding_sites.csv")


def _load_runner(engine_root: Path):
    """Import only the pinned stdlib-only runner, without changing sys.path."""
    path = engine_root / "integrations/mcp/runner.py"
    spec = importlib.util.spec_from_file_location("zincsight_batch_pinned_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _validate_inputs(structures, histidine_rotamers, minimum_probability,
                     max_sites, cpu_workers, timeout):
    if not isinstance(structures, list) or not 1 <= len(structures) <= MAX_STRUCTURES:
        raise ValueError(f"Provide between 1 and {MAX_STRUCTURES} structures.")
    if not isinstance(histidine_rotamers, bool):
        raise ValueError("histidine_rotamers must be a boolean.")
    if (isinstance(minimum_probability, bool)
            or not isinstance(minimum_probability, (int, float))
            or not math.isfinite(minimum_probability)
            or not 0 <= minimum_probability <= 1):
        raise ValueError("minimum_probability must be between 0 and 1.")
    for name, value, upper in (("max_sites", max_sites, 100), ("cpu_workers", cpu_workers, 8)):
        if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= upper:
            raise ValueError(f"{name} must be an integer between 1 and {upper}.")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be a positive number of seconds.")
    seen = set()
    total = 0
    validated = []
    for item in structures:
        if not isinstance(item, dict):
            raise ValueError("Every structure must have a key, source, structure_text and file_format.")
        key, source = item.get("key"), item.get("source")
        content, file_format = item.get("structure_text"), item.get("file_format")
        if not isinstance(key, str) or not KEY_PATTERN.fullmatch(key) or key in seen:
            raise ValueError("Each structure key must be unique and contain 1–64 letters, digits, underscores or hyphens, starting with a letter or digit.")
        seen.add(key)
        if not isinstance(source, str) or not source.strip():
            raise ValueError(f"Structure {key}: source must be a nonempty string.")
        if file_format not in ("pdb", "cif", "mmcif"):
            raise ValueError(f"Structure {key}: file_format must be pdb, cif or mmcif.")
        file_format = "cif" if file_format == "mmcif" else file_format
        if not isinstance(content, str) or not content.strip() or "\x00" in content:
            raise ValueError(f"Structure {key}: provide nonempty text PDB/mmCIF content.")
        size = len(content.encode("utf-8"))
        if size > MAX_BYTES:
            raise ValueError(f"Structure {key} exceeds the 5 MB limit.")
        total += size
        if total > MAX_TOTAL_BYTES:
            raise ValueError("Combined structures exceed the 100 MB limit.")
        if file_format == "pdb" and not any(line.startswith("ATOM  ") for line in content.splitlines()):
            raise ValueError(f"Structure {key}: PDB input must contain protein ATOM records.")
        if file_format == "cif" and "_atom_site." not in content:
            raise ValueError(f"Structure {key}: mmCIF input must contain an atom_site table.")
        validated.append({"key": key, "source": source, "structure_text": content, "file_format": file_format})
    return validated


def _parse_batch_sites(table: Path, keys: list[str], minimum_probability: float,
                       parse_sites) -> dict[str, list[dict]]:
    """Split native structure IDs, then use the pinned parser for every score."""
    grouped = {key: [] for key in keys}
    with table.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        if not fields or "structure_id" not in fields:
            raise ValueError("ZincSight batch output is missing structure_id.")
        for row in reader:
            key = row["structure_id"]
            if key not in grouped:
                raise ValueError(f"ZincSight batch output contains an unknown structure_id: {key!r}.")
            grouped[key].append(row)
    # All files live inside the disposable worker directory. Header-only groups
    # represent inputs with no sites, including mixed positive/negative batches.
    for key, rows in grouped.items():
        path = table.parent / f"parsed-{key}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        grouped[key] = parse_sites(path, minimum_probability)
    return grouped


def _emit_diagnostics(work: Path):
    for name in ("worker.log", "database.log"):
        path = work / name
        if path.is_file():
            with path.open("rb") as stream:
                stream.seek(max(0, path.stat().st_size - 32768))
                tail = stream.read().decode("utf-8", errors="replace")
            print(f"--- ZincSight batch {name} (last 32 KiB) ---\n{tail}", file=sys.stderr)


def run_predictions(structures: list[dict], *, histidine_rotamers: bool,
                    minimum_probability: float, max_sites: int, cpu_workers: int,
                    engine_root: Path = Path("/opt/zincsight"), timeout: int = 1200) -> list[dict]:
    """Predict all inputs once and return one result per input, in input order.

    Inputs carry a unique safe ``key``, display ``source``, ``structure_text`` and
    ``file_format``. Thresholding and truncation are independent for each input.
    The original first model is normalized by the pinned scientific worker.
    """
    structures = _validate_inputs(structures, histidine_rotamers, minimum_probability,
                                 max_sites, cpu_workers, timeout)
    engine_root = Path(engine_root).resolve()
    runner = _load_runner(engine_root)
    with tempfile.TemporaryDirectory(prefix="zincsight-batch-") as directory:
        work = Path(directory)
        try:
            shutil.copytree(engine_root / "src", work / "src",
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            source = engine_root / "integrations/mcp"
            shutil.copy2(source / "worker.py", work / "pinned_worker.py")
            shutil.copy2(source / "embedded_database.py", work / "embedded_database.py")
            shutil.copy2(Path(__file__).with_name("batch_worker.py"), work / "batch_worker.py")
            inputs = work / "inputs"
            inputs.mkdir()
            manifest = []
            for item in structures:
                filename = f"inputs/{item['key']}.{item['file_format']}"
                (work / filename).write_text(item["structure_text"], encoding="utf-8")
                manifest.append({"key": item["key"], "source": item["source"],
                                 "file_format": item["file_format"], "path": filename})
            (work / "batch-inputs.json").write_text(json.dumps(manifest), encoding="utf-8")
            environment = os.environ.copy()
            environment.update(
                PYTHON_DOTENV_DISABLED="1", DEBUGGING_TEMP="false", KEEP_TEMP_TABLES_TEMP="false",
                FILTER_PROB_AFTER_COMPRESSION="false", EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD="0.0",
                OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg",
                ZINCSIGHT_DATABASE_BACKEND="pglite", ZINCSIGHT_MAX_RESIDUES=str(MAX_RESIDUES),
                ZINCSIGHT_PGLITE_LAUNCHER=str(source / "pglite/server.mjs"))
            with (work / "worker.log").open("w", encoding="utf-8") as log:
                process = subprocess.Popen(
                    [sys.executable, str(work / "batch_worker.py"), str(int(histidine_rotamers)), str(cpu_workers)],
                    cwd=work, env=environment, stdout=log, stderr=log, start_new_session=True)
                try:
                    returncode = process.wait(timeout=timeout)
                except subprocess.TimeoutExpired as exc:
                    raise RuntimeError(f"ZincSight batch exceeded the {timeout}-second limit; use fewer or smaller structures.") from exc
                finally:
                    # Kill the entire session, including any database or pool
                    # children that survive a failure, before temporary cleanup.
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
            status_path = work / "status.json"
            if returncode != 0 or not status_path.is_file():
                raise RuntimeError(f"ZincSight batch worker failed (exit {returncode}); see worker/database logs for the affected structure and cause.")
            status = json.loads(status_path.read_text(encoding="utf-8"))
            counts = status.get("protein_residue_counts", {})
            keys = [item["key"] for item in structures]
            if (set(counts) != set(keys) or not isinstance(status.get("has_sites"), bool)
                    or any(isinstance(counts[key], bool) or not isinstance(counts[key], int)
                           or not 3 <= counts[key] <= MAX_RESIDUES for key in keys)):
                raise RuntimeError("ZincSight batch worker returned invalid per-structure validation metadata.")
            sites = (_parse_batch_sites(work / RESULT_TABLE, keys, minimum_probability, runner.parse_sites)
                     if status["has_sites"] else {key: [] for key in keys})
            return [{
                "method": "ZincSight", "key": item["key"], "source": item["source"],
                "histidine_rotamers": histidine_rotamers, "minimum_probability": minimum_probability,
                "site_count": len(sites[item["key"]]), "sites": sites[item["key"]][:max_sites],
                "truncated": len(sites[item["key"]]) > max_sites,
                "model_policy": MODEL_POLICY, "database_backend": "pglite",
                "protein_residue_count": counts[item["key"]],
                "limitations": runner.LIMITATIONS, "software_url": "https://github.com/MECHTI1/ZincSight",
            } for item in structures]
        except Exception:
            _emit_diagnostics(work)
            raise
