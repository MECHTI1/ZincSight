"""Exercise the mounted-source image offline and retain exact scientific parity checks."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PARITY_KEYS = ("sites", "site_count", "truncated", "histidine_rotamers",
               "minimum_probability", "model_policy")

# The host need not have BioPython or PyMOL; this executes in the image.
ROUNDTRIP_CHECK = r'''
import csv, json
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pymol2

out = Path('/app/out')
batch = json.loads((out / 'prediction.json').read_text())
with (out / 'structures.csv').open() as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == batch['structure_count'] == len(batch['structures'])
for row, result in zip(rows, batch['structures']):
    assert row['structure_id'] == result['key']
    assert int(row['returned_site_count']) == len(result['sites'])
    for field in ('structure_file', 'mmcif_file', 'pymol_session', 'prediction_json', 'site_annotations'):
        path = (out / row[field]).resolve()
        assert path.is_relative_to(out.resolve()) and path.is_file() and path.stat().st_size
    annotations = json.loads((out / row['site_annotations']).read_text())
    assert annotations['sites'] == result['sites']
    assert annotations['site_count'] == len(result['sites'])
    assert '10.1002/pro.70350' in annotations['citation']
    cif = MMCIF2Dict(str(out / row['mmcif_file']))
    chains = cif['_atom_site.auth_asym_id']
    marker_chain = annotations['predicted_marker_chain']
    marker_indexes = [i for i, chain in enumerate(chains) if chain == marker_chain]
    assert len(marker_indexes) == len(result['sites'])
    original_atom_count = len(cif['_atom_site.id']) - len(marker_indexes)
    assert len(set(cif['_atom_site.pdbx_PDB_model_num'])) == 1
    for index, site in zip(marker_indexes, result['sites']):
        assert cif['_atom_site.group_PDB'][index] == 'HETATM'
        assert cif['_atom_site.type_symbol'][index].upper() == 'ZN'
        assert [float(cif['_atom_site.Cartn_' + axis][index]) for axis in 'xyz'] == site['coordinates_angstrom']
    pdb = PDBParser(QUIET=True).get_structure('export', str(out / result['artifacts']['pdb']))
    assert len(list(pdb.get_models())) == 1
    assert len(list(pdb.get_atoms())) == original_atom_count + len(result['sites'])
    if result['sites']:
        markers = list(pdb[0][annotations['pdb_predicted_marker_chain']].get_atoms())
        assert len(markers) == len(result['sites'])
        for atom, site in zip(markers, result['sites']):
            assert all(abs(float(a) - b) <= 0.00051 for a, b in zip(atom.coord, site['coordinates_angstrom']))
    with pymol2.PyMOL() as session:
        cmd = session.cmd
        cmd.load(str(out / row['pymol_session']))
        assert cmd.count_states('original_structure') == 1
        assert cmd.count_atoms('original_structure') == original_atom_count
        original_metals = sum(symbol.upper() in ('ZN', 'FE', 'MN', 'CO', 'CU', 'NI', 'MG', 'CA')
                              for i, symbol in enumerate(cif['_atom_site.type_symbol'])
                              if i not in marker_indexes)
        assert cmd.count_atoms('experimental_metals') == original_metals
        if result['sites']:
            assert cmd.count_atoms('predicted_sites') == len(result['sites'])
            assert cmd.count_atoms('coordinating_residues') > 0
            for atom, site in zip(cmd.get_model('predicted_sites').atom, result['sites']):
                assert all(abs(a - b) < 0.0001 for a, b in zip(atom.coord, site['coordinates_angstrom']))
            labels = []
            cmd.iterate('predicted_sites', 'labels.append(label)', space={'labels': labels})
            assert all('Predicted site' in label and 'p=' in label for label in labels)
        else:
            assert 'predicted_sites' not in cmd.get_names('objects')
    print('PASS structure/PyMOL roundtrip:', result['key'], len(result['sites']), 'markers', flush=True)
'''


def run_checked(command, log, *, timeout=360):
    process = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    log.write_text(process.stdout + process.stderr)
    if process.returncode:
        raise RuntimeError(f"Command failed; see {log}:\n{process.stdout}\n{process.stderr}")
    return process


def check_csv(rows, results):
    expected = [(result, site) for result in results for site in result["sites"]]
    assert len(rows) == len(expected)
    for row, (result, site) in zip(rows, expected):
        assert row["structure_id"] == result["key"]
        assert float(row["probability"]) == site["probability"]
        assert [float(row[k]) for k in ("x_angstrom", "y_angstrom", "z_angstrom")] == site["coordinates_angstrom"]
        assert json.loads(row["coordinating_residues"]) == site["coordinating_residues"]
        assert row["template_id"] == str(site["template_id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", default="zincsight-tamarind:test")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="zincsight-tamarind-test-") as temporary:
        app = Path(temporary)
        (app / "inputs").mkdir()
        shutil.copy2(ROOT / "run.sh", app / "run.sh")
        # Tamarind replaces /app: mount every adapter module, keeping /opt usable.
        shutil.copytree(ROOT / "integrations/tamarind", app / "integrations/tamarind",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name in ("1KLS.pdb", "1KLS.cif", "2A0S-assembly1.cif"):
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{name}", timeout=60) as response:
                raw = response.read(5_000_001)
            assert len(raw) <= 5_000_000
            (app / "inputs" / name).write_bytes(raw)
        (app / "inputs/invalid.pdb").write_text("This is not a protein structure.\n")
        (app / "inputs/empty_sites.pdb").write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
            "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 90.00           C  \n"
            "ATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00 90.00           C  \nEND\n")
        (app / "inputs/duplicate").mkdir()
        shutil.copy2(app / "inputs/1KLS.pdb", app / "inputs/duplicate/1KLS.pdb")
        hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in (app / "inputs").glob("*") if p.is_file()}
        (args.output / "fixture_sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")
        spec = importlib.util.spec_from_file_location("prepare_inputs", ROOT / "integrations/tamarind/prepare_inputs.py")
        preparer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(preparer)
        manifest = preparer.prepare_inputs(app / "inputs/prepared.zip",
                                           sources=[app / "inputs/1KLS.cif", app / "inputs/empty_sites.pdb"])
        base = ["docker", "run", "--rm", "--network", "none", "--cpus", "4", "--memory", "12g",
                "--user", f"{os.getuid()}:{os.getgid()}", "--volume", f"{app}:/app"]
        resources = ["-e", "ZINCSIGHT_CPU_BUDGET=4", "-e", "ZINCSIGHT_MEMORY_GIB=12", "-e", "cpu_workers=4"]
        references = {}

        def reference(name, rotamers, limit):
            key = (name, rotamers, limit)
            if key not in references:
                command = ["python", "-m", "integrations.mcp.predict", "--database", "pglite",
                           "--structure", f"/app/inputs/{name}", "--format", "pdb" if name.endswith(".pdb") else "cif",
                           "--max-sites", limit, "--minimum-probability", "0.0001"]
                if rotamers == "true":
                    command.append("--histidine-rotamers")
                log = args.output / (name.replace("/", "_") + f"-{rotamers}-{limit}-reference.log")
                direct = run_checked(base + ["--workdir", "/opt/zincsight", args.image] + command, log)
                references[key] = json.loads(direct.stdout)
            return references[key]

        def case(label, names, environment, rotamers="true", limit="100"):
            output = args.output / label
            output.mkdir(exist_ok=True)
            shutil.rmtree(app / "out", ignore_errors=True)
            settings = ["-e", f"histidine_rotamers={rotamers}", "-e", f"max_sites={limit}"]
            for key, value in environment.items():
                settings.extend(["-e", f"{key}={value}"])
            run_checked(base + resources + settings + [args.image], output / "container.log")
            batch = json.loads((app / "out/prediction.json").read_text())
            results = batch["structures"]
            assert batch["structure_count"] == len(results) == len(names)
            assert len({result["key"] for result in results}) == len(names)
            assert batch["performance"]["database_instances_per_batch"] == 1
            for name, result in zip(names, results):
                baseline = reference("1KLS.pdb" if name == "duplicate/1KLS.pdb" else name, rotamers, limit)
                for key in PARITY_KEYS:
                    assert result[key] == baseline[key], (label, name, key)
                assert result["input_sha256"] == hashlib.sha256((app / "inputs" / name).read_bytes()).hexdigest()
                assert json.loads((app / "out" / result["prediction_json"]).read_text()) == result
            with (app / "out/sites.csv").open() as stream:
                check_csv(list(csv.DictReader(stream)), results)
            assert batch["returned_site_count"] == sum(len(result["sites"]) for result in results)
            assert batch["site_count"] == sum(result["site_count"] for result in results)
            for filename in ("CITATION.txt", "CITATION.bib", "citation.json"):
                assert "10.1002/pro.70350" in (app / "out" / filename).read_text()
            with zipfile.ZipFile(app / "out/zincsight_results.zip") as archive:
                expected = {p.relative_to(app / "out").as_posix() for p in (app / "out").rglob("*")
                            if p.is_file() and p.name != "zincsight_results.zip"}
                assert set(archive.namelist()) == expected
                assert all(archive.read(name) == (app / "out" / name).read_bytes() for name in expected)
            run_checked(base + ["--entrypoint", "python", args.image, "-c", ROUNDTRIP_CHECK], output / "roundtrip.log")
            shutil.copytree(app / "out", output / "out", dirs_exist_ok=True)
            print(f"PASS {label}: {len(results)} structures, exact wrapper agreement and session roundtrip", flush=True)
            return batch

        for label, name, rotamers, limit in [
                ("pdb-rotamers", "1KLS.pdb", "true", "100"),
                ("cif-rotamers", "1KLS.cif", "true", "100"),
                ("cif-no-rotamers", "1KLS.cif", "false", "100"),
                ("truncated", "2A0S-assembly1.cif", "true", "1"),
                ("no-sites", "empty_sites.pdb", "true", "100")]:
            batch = case(label, [name], {"structure": f"/app/inputs/{name}"}, rotamers, limit)
            result = batch["structures"][0]
            # Preserve the flat single-input JSON compatibility and exact baseline checks.
            assert all(batch[key] == result[key] for key in PARITY_KEYS)
            if label == "truncated":
                assert result["truncated"] and result["site_count"] > 1 and len(result["sites"]) == 1
            elif label == "no-sites":
                assert result["site_count"] == 0 and not result["truncated"]
            else:
                assert result["site_count"] >= 1

        names = ["1KLS.pdb", "duplicate/1KLS.pdb", "empty_sites.pdb"]
        batch = case("batch-list", names, {"structures": json.dumps([f"/app/inputs/{name}" for name in names])})
        assert batch["structures"][0]["source"] == batch["structures"][1]["source"] == "1KLS.pdb"
        assert batch["structures"][2]["sites"] == []
        batch = case("batch-prepared-archive", ["1KLS.cif", "empty_sites.pdb"],
                     {"batch_archive": "/app/inputs/prepared.zip"})
        for result, entry in zip(batch["structures"], manifest["structures"]):
            assert result["input_provenance"] == entry

        shutil.rmtree(app / "out", ignore_errors=True)
        invalid = subprocess.run(base + resources + ["-e", "structure=/app/inputs/invalid.pdb", args.image],
                                 text=True, capture_output=True, timeout=60)
        (args.output / "invalid.log").write_text(invalid.stdout + invalid.stderr)
        assert invalid.returncode != 0 and not (app / "out/prediction.json").exists()
        print("PASS invalid structure: nonzero exit, no prediction JSON", flush=True)


if __name__ == "__main__":
    main()

