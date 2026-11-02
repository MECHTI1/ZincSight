"""Run the real mounted-source container with networking disabled (Docker required)."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", default="zincsight-tamarind:test")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="zincsight-tamarind-test-") as temporary:
        app = Path(temporary)
        (app / "inputs").mkdir()
        (app / "integrations/tamarind").mkdir(parents=True)
        shutil.copy2(ROOT / "run.sh", app / "run.sh")
        shutil.copy2(ROOT / "integrations/tamarind/predict.py", app / "integrations/tamarind/predict.py")
        # Only the adapter is mounted: the scientific engine and all dependencies
        # must be available from /opt even when /app is replaced by Tamarind.
        for name in ("1KLS.pdb", "1KLS.cif", "2A0S-assembly1.cif"):
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{name}", timeout=60) as response:
                (app / "inputs" / name).write_bytes(response.read())
        (app / "inputs/invalid.pdb").write_text("This is not a protein structure.\n")
        (app / "inputs/empty_sites.pdb").write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
            "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 90.00           C  \n"
            "ATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00 90.00           C  \nEND\n")
        base = ["docker", "run", "--rm", "--network", "none", "--cpus", "2", "--memory", "8g",
                "--user", f"{os.getuid()}:{os.getgid()}", "--volume", f"{app}:/app"]
        cases = [("pdb-rotamers", "1KLS.pdb", "true", "100"),
                 ("cif-rotamers", "1KLS.cif", "true", "100"),
                 ("cif-no-rotamers", "1KLS.cif", "false", "100"),
                 ("truncated", "2A0S-assembly1.cif", "true", "1"),
                 ("no-sites", "empty_sites.pdb", "true", "100")]
        for label, name, rotamers, limit in cases:
            output = args.output / label
            output.mkdir()
            shutil.rmtree(app / "out", ignore_errors=True)
            environment = ["-e", f"structure=/app/inputs/{name}", "-e", f"histidine_rotamers={rotamers}",
                           "-e", f"max_sites={limit}"]
            process = subprocess.run(base + environment + [args.image], text=True, capture_output=True, timeout=240)
            (output / "container.log").write_text(process.stdout + process.stderr)
            if process.returncode:
                raise RuntimeError(f"Container failed for {label}: {process.stdout}\n{process.stderr}")
            result = json.loads((app / "out/prediction.json").read_text())
            shutil.copytree(app / "out", output / "out")
            # The same exact input and settings must match the established wrapper.
            command = ["python", "-m", "integrations.mcp.predict", "--database", "pglite",
                       "--structure", f"/app/inputs/{name}", "--format", "pdb" if name.endswith(".pdb") else "cif",
                       "--max-sites", limit, "--minimum-probability", "0.0001"]
            if rotamers == "true":
                command.append("--histidine-rotamers")
            direct = subprocess.run(base + ["--workdir", "/opt/zincsight", args.image] + command,
                                    text=True, capture_output=True, check=True, timeout=240)
            reference = json.loads(direct.stdout)
            for key in ("sites", "site_count", "truncated", "histidine_rotamers", "minimum_probability", "model_policy"):
                assert result[key] == reference[key], (label, key)
            assert result["input_sha256"] == hashlib.sha256((app / "inputs" / name).read_bytes()).hexdigest()
            with (app / "out/sites.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            assert len(rows) == result["returned_site_count"] == len(result["sites"])
            for row, site in zip(rows, result["sites"]):
                assert float(row["probability"]) == site["probability"]
                assert [float(row[k]) for k in ("x_angstrom", "y_angstrom", "z_angstrom")] == site["coordinates_angstrom"]
                assert json.loads(row["coordinating_residues"]) == site["coordinating_residues"]
            if label == "truncated":
                assert result["truncated"] and result["site_count"] > 1 and len(rows) == 1
            elif label == "no-sites":
                assert result["site_count"] == 0 and not result["truncated"]
            else:
                assert result["site_count"] >= 1
            print(f"PASS {label}: {len(rows)}/{result['site_count']} sites, exact wrapper agreement", flush=True)
        shutil.rmtree(app / "out", ignore_errors=True)
        invalid = subprocess.run(base + ["-e", "structure=/app/inputs/invalid.pdb", args.image],
                                 text=True, capture_output=True, timeout=60)
        assert invalid.returncode != 0 and not (app / "out/prediction.json").exists()
        print("PASS invalid structure: nonzero exit, no prediction JSON", flush=True)


if __name__ == "__main__":
    main()
