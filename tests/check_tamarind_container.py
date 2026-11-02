"""Exercise the mounted-source image offline and retain exact scientific parity checks."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import tempfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PARITY_KEYS = ("sites", "site_count", "truncated", "histidine_rotamers",
               "minimum_probability", "model_policy")
PLATFORM_LOG_PREFIX = "Platform output log opened before the adapter.\n"
PLATFORM_STATE = '{"job_state":"running","owner":"platform"}\n'
PLATFORM_WRAPPER = r'''set -e
trap 'status=$?; printf "%s\n" "$status" > /app/out/.exit_code' EXIT
bash /app/run.sh >> /app/out/output.log 2>&1
'''

# The host need not have BioPython or PyMOL; this executes in the image.
ROUNDTRIP_CHECK = r'''
from collections import Counter, defaultdict
import csv, hashlib, json
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pymol2


def unpacked_atoms(model):
    # get_atoms() selects one conformer per DisorderedAtom/DisorderedResidue.
    # PDBIO exports all conformers, so inspect the same complete atom set here.
    return [atom for chain in model.get_chains()
            for residue in chain.get_unpacked_list()
            for atom in residue.get_unpacked_list()]


def clean(value):
    return '' if value in ('.', '?', ' ') else str(value)


def cif_records(data, indexes, annotations=None):
    records = defaultdict(list)
    for index in indexes:
        chain = clean(data['_atom_site.auth_asym_id'][index])
        number = str(int(data['_atom_site.auth_seq_id'][index]))
        insertion = clean(data['_atom_site.pdbx_PDB_ins_code'][index])
        resi = number + insertion
        if annotations is not None:
            mapped = annotations['pdb_residue_mapping'].get(f'{chain}_{resi}')
            if mapped:
                chain, resi = clean(mapped['chain']), mapped['resi']
            else:
                chain = clean(annotations['pdb_chain_mapping'][chain or ' '])
        key = (chain, resi, data['_atom_site.auth_comp_id'][index].strip(),
               data['_atom_site.auth_atom_id'][index].strip(),
               clean(data['_atom_site.label_alt_id'][index]),
               data['_atom_site.type_symbol'][index].upper(),
               data['_atom_site.group_PDB'][index] == 'HETATM')
        records[key].append(tuple(float(data['_atom_site.Cartn_' + axis][index]) for axis in 'xyz'))
    return records


def pdb_records(atoms):
    records = defaultdict(list)
    for atom in atoms:
        residue = atom.get_parent()
        key = (clean(residue.get_parent().id), str(residue.id[1]) + clean(residue.id[2]),
               residue.resname.strip(), atom.name.strip(), clean(atom.altloc), atom.element.upper(),
               residue.id[0] != ' ')
        records[key].append(tuple(float(value) for value in atom.coord))
    return records


def pymol_records(atoms):
    records = defaultdict(list)
    for atom in atoms:
        key = (clean(atom.chain), atom.resi, atom.resn.strip(), atom.name.strip(), clean(atom.alt),
               atom.symbol.upper(), bool(atom.hetatm))
        records[key].append(tuple(atom.coord))
    return records


def assert_records(actual, expected, tolerance):
    # Matching identities prevents atom losses/replacements or altloc/chain swaps
    # from hiding behind an equal total atom count.
    assert Counter({key: len(value) for key, value in actual.items()}) == Counter(
        {key: len(value) for key, value in expected.items()}), 'Atom identifiers/alternate locations differ'
    for key in expected:
        for observed, reference in zip(sorted(actual[key]), sorted(expected[key])):
            assert all(abs(a - b) <= tolerance for a, b in zip(observed, reference)), (key, observed, reference)


out = Path('/app/details_bundle')
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
    marker_index_set = set(marker_indexes)
    original_indexes = [i for i in range(len(cif['_atom_site.id'])) if i not in marker_index_set]
    original_atom_count = len(original_indexes)
    assert len(set(cif['_atom_site.pdbx_PDB_model_num'])) == 1
    original_records = cif_records(cif, original_indexes)
    # These fixture names are also staged directly, including prepared ZIP inputs.
    source = Path('/app/inputs') / result['source']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == result['input_sha256']
    if source.suffix.lower() in ('.cif', '.mmcif'):
        source_cif = MMCIF2Dict(str(source))
        models = source_cif.get('_atom_site.pdbx_PDB_model_num', ['1'] * len(source_cif['_atom_site.id']))
        first_model_indexes = [i for i, model in enumerate(models) if model == models[0]]
        # Compare original atom identities and coordinate strings, not unrelated
        # structure metadata copied alongside the atom table.
        identity_columns = ('group_PDB', 'type_symbol', 'label_atom_id', 'auth_atom_id',
                            'label_comp_id', 'auth_comp_id', 'label_asym_id', 'auth_asym_id',
                            'label_seq_id', 'auth_seq_id', 'label_alt_id', 'pdbx_PDB_ins_code',
                            'Cartn_x', 'Cartn_y', 'Cartn_z')
        for field in identity_columns:
            column = '_atom_site.' + field
            if column in source_cif:
                values = source_cif[column]
                assert [cif[column][i] for i in original_indexes] == [values[i] for i in first_model_indexes], column
    else:
        source_model = next(PDBParser(QUIET=True).get_structure('source', str(source)).get_models())
        assert_records(original_records, pdb_records(unpacked_atoms(source_model)), 0.00001)
    for index, site in zip(marker_indexes, result['sites']):
        assert cif['_atom_site.group_PDB'][index] == 'HETATM'
        assert cif['_atom_site.type_symbol'][index].upper() == 'ZN'
        assert [float(cif['_atom_site.Cartn_' + axis][index]) for axis in 'xyz'] == site['coordinates_angstrom']
    pdb = PDBParser(QUIET=True).get_structure('export', str(out / result['artifacts']['pdb']))
    assert len(list(pdb.get_models())) == 1
    pdb_atoms = unpacked_atoms(pdb[0])
    assert len(pdb_atoms) == original_atom_count + len(result['sites'])
    pdb_original_atoms = [atom for atom in pdb_atoms
                          if atom.get_parent().get_parent().id != annotations['pdb_predicted_marker_chain']]
    assert_records(pdb_records(pdb_original_atoms), cif_records(cif, original_indexes, annotations), 0.00051)
    if result['sites']:
        markers = [atom for atom in pdb_atoms
                   if atom.get_parent().get_parent().id == annotations['pdb_predicted_marker_chain']]
        assert len(markers) == len(result['sites'])
        for atom, site in zip(markers, result['sites']):
            assert all(abs(float(a) - b) <= 0.00051 for a, b in zip(atom.coord, site['coordinates_angstrom']))
    with pymol2.PyMOL() as session:
        cmd = session.cmd
        cmd.load(str(out / row['pymol_session']))
        assert cmd.count_states('original_structure') == 1
        # mmCIF polymer metadata may add coordinate-less sequence placeholders;
        # count_atoms(state=0) includes those, unlike the exported coordinate model.
        assert cmd.count_atoms('original_structure', state=1) == original_atom_count
        assert_records(pymol_records(cmd.get_model('original_structure').atom), original_records, 0.0001)
        original_metals = sum(symbol.upper() in ('ZN', 'FE', 'MN', 'CO', 'CU', 'NI', 'MG', 'CA')
                              for i, symbol in enumerate(cif['_atom_site.type_symbol'])
                              if i not in marker_indexes)
        assert cmd.count_atoms('experimental_metals', state=1) == original_metals
        if result['sites']:
            assert cmd.count_atoms('predicted_sites', state=1) == len(result['sites'])
            assert cmd.count_atoms('coordinating_residues', state=1) > 0
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


def prepare_platform_outputs(output):
    """Create runner-owned output before Docker and emulate service-owned files."""
    shutil.rmtree(output, ignore_errors=True)
    output.mkdir()
    (output / "output.log").write_text(PLATFORM_LOG_PREFIX)
    (output / ".exit_code").write_text("pending\n")
    (output / ".platform_state").mkdir()
    (output / ".platform_state/state.json").write_text(PLATFORM_STATE)
    return {name: (path.stat().st_dev, path.stat().st_ino)
            for name in ("output.log", ".exit_code", ".platform_state", ".platform_state/state.json")
            for path in (output / name,)}


def check_platform_outputs(output, canaries, *, exit_code=0):
    for name, identity in canaries.items():
        info = (output / name).stat()
        assert (info.st_dev, info.st_ino) == identity, f"Platform path was moved or replaced: {name}"
    assert (output / ".platform_state/state.json").read_text() == PLATFORM_STATE
    assert (output / ".exit_code").read_text() == f"{exit_code}\n"
    log = (output / "output.log").read_text()
    assert log.startswith(PLATFORM_LOG_PREFIX), "The open service log was lost or truncated"
    if exit_code == 0:
        assert "Processed " in log and "PyMOL sessions are ready." in log, "Final adapter output missed the open service log"


def check_public_outputs(output, batch, details, *, service_canaries=None):
    """Check the uncluttered surface, then recover the complete detailed layout."""
    assert batch["integration_version"] == "0.2.3"
    assert batch["output_layout"] == "summary_and_bundle"
    assert batch["result_bundle"] == "ZincSight_results.zip"
    with (output / "00_Summary.csv").open() as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ["structure", "sites", "top_probability", "status", "structure_file"]
        summary = list(reader)
    with (output / "01_Predicted_sites.csv").open() as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ["structure", "site", "probability", "coordinating_residues",
                                     "template_id", "template_rmsd_angstrom"]
        site_rows = list(reader)
    assert len(summary) == batch["structure_count"] == len(batch["structures"])
    expected_files = {"00_Summary.csv", "01_Predicted_sites.csv", "prediction.json",
                      "CITATION.txt", "ZincSight_results.zip"}
    expected_sites = []
    for row, result in zip(summary, batch["structures"]):
        assert row["structure"] == result["display_name"]
        assert int(row["sites"]) == result["site_count"]
        expected_status = "No sites predicted" if not result["sites"] else "Truncated" if result["truncated"] else "Complete"
        assert row["status"] == expected_status
        if result["sites"]:
            assert float(row["top_probability"]) == max(site["probability"] for site in result["sites"])
        else:
            assert row["top_probability"] == ""
        viewer = row["structure_file"]
        assert Path(viewer).name == viewer and re.fullmatch(r"\d{2,3}_.+_ZincSight\.(pdb|cif)", viewer)
        session = Path(viewer).with_suffix(".pse").name
        assert not {viewer, session} & expected_files, "Public input filenames collide"
        expected_files.update((viewer, session))
        artifacts = result["artifacts"]
        assert artifacts["structure"] == viewer
        assert artifacts["pymol_session"] == session
        assert artifacts["details_bundle"] == "ZincSight_results.zip"
        assert artifacts["pdb" if viewer.endswith(".pdb") else "mmcif"] == viewer
        assert set(artifacts.values()) == {viewer, session, "ZincSight_results.zip"}
        for rank, site in enumerate(result["sites"], 1):
            expected_sites.append((row["structure"], rank, site))
    assert len(site_rows) == len(expected_sites)
    for row, (structure, rank, site) in zip(site_rows, expected_sites):
        assert row["structure"] == structure and int(row["site"]) == rank
        assert float(row["probability"]) == site["probability"]
        assert row["coordinating_residues"] == ", ".join(site["coordinating_residues"])
        assert row["template_id"] == str(site["template_id"])
        assert float(row["template_rmsd_angstrom"]) == site["template_rmsd_angstrom"]
    service_roots = {name.split("/", 1)[0] for name in service_canaries or ()}
    assert {path.name for path in output.iterdir()} == expected_files | service_roots
    assert all((output / name).is_file() and (output / name).stat().st_size for name in expected_files)
    if service_canaries:
        check_platform_outputs(output, service_canaries)
    shutil.rmtree(details, ignore_errors=True)
    details.mkdir()
    with zipfile.ZipFile(output / "ZincSight_results.zip") as archive:
        members = archive.infolist()
        names = {info.filename for info in members}
        assert len(names) == len(members), "Duplicate bundle member names"
        assert not {"zincsight_results.zip", "ZincSight_results.zip"} & names, "Recursive results bundle"
        assert not any(name == root or name.startswith(root + "/")
                       for name in names for root in service_roots), "Platform files entered the results bundle"
        assert {"structures.csv", "sites.csv", "prediction.json", "CITATION.txt", "CITATION.bib",
                "citation.json", "00_Summary.csv", "01_Predicted_sites.csv", "public_manifest.json"} <= names
        for info in members:
            destination = (details / info.filename).resolve()
            assert destination.is_relative_to(details.resolve()) and not stat.S_ISLNK(info.external_attr >> 16)
            if info.is_dir():
                destination.mkdir(parents=True, exist_ok=True)
            else:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(archive.read(info))  # Also verifies each member's CRC.
        for filename in ("00_Summary.csv", "01_Predicted_sites.csv", "CITATION.txt"):
            assert archive.read(filename) == (output / filename).read_bytes()
        bundled = json.loads(archive.read("prediction.json"))
        assert len(bundled["structures"]) == len(batch["structures"])
        for row, result, original in zip(summary, batch["structures"], bundled["structures"]):
            for key in PARITY_KEYS + ("key", "source", "input_sha256", "input_provenance"):
                assert result[key] == original[key], key
            assert result["prediction_json_in_bundle"] == original["prediction_json"]
            assert json.loads(archive.read(original["prediction_json"])) == original
            for name in result["bundle_members"].values():
                assert name in names
            for kind in ("pdb" if row["structure_file"].endswith(".pdb") else "mmcif", "pymol_session"):
                assert (output / result["artifacts"][kind]).read_bytes() == archive.read(original["artifacts"][kind])
    return bundled


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
                "--user", f"{os.getuid()}:{os.getgid()}", "--volume", f"{app}:/app",
                "--volume", f"{app / 'out'}:/app/out"]
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
            # A separate output mount cannot be renamed/replaced. Create it as
            # the runner first; Docker would otherwise create a root-owned path.
            canaries = prepare_platform_outputs(app / "out")
            settings = ["-e", f"histidine_rotamers={rotamers}", "-e", f"max_sites={limit}"]
            for key, value in environment.items():
                settings.extend(["-e", f"{key}={value}"])
            try:
                run_checked(base + resources + settings + ["--entrypoint", "bash", args.image,
                                                          "-c", PLATFORM_WRAPPER], output / "container.log")
            finally:
                # Preserve the service log even when the adapter/container fails.
                shutil.copytree(app / "out", output / "out", dirs_exist_ok=True)
            batch = json.loads((app / "out/prediction.json").read_text())
            results = batch["structures"]
            assert batch["structure_count"] == len(results) == len(names)
            assert len({result["key"] for result in results}) == len(names)
            assert batch["performance"]["database_instances_per_batch"] == 1
            bundled = check_public_outputs(app / "out", batch, app / "details_bundle", service_canaries=canaries)
            for name, result in zip(names, results):
                baseline = reference("1KLS.pdb" if name == "duplicate/1KLS.pdb" else name, rotamers, limit)
                for key in PARITY_KEYS:
                    assert result[key] == baseline[key], (label, name, key)
                assert result["input_sha256"] == hashlib.sha256((app / "inputs" / name).read_bytes()).hexdigest()
            with (app / "details_bundle/sites.csv").open() as stream:
                check_csv(list(csv.DictReader(stream)), bundled["structures"])
            assert batch["returned_site_count"] == sum(len(result["sites"]) for result in results)
            assert batch["site_count"] == sum(result["site_count"] for result in results)
            for filename in ("CITATION.txt", "CITATION.bib", "citation.json"):
                assert "10.1002/pro.70350" in (app / "details_bundle" / filename).read_text()
            run_checked(base + ["--entrypoint", "python", args.image, "-c", ROUNDTRIP_CHECK], output / "roundtrip.log")
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

        canaries = prepare_platform_outputs(app / "out")
        invalid = subprocess.run(base + resources + ["-e", "structure=/app/inputs/invalid.pdb",
                                                    "--entrypoint", "bash", args.image, "-c", PLATFORM_WRAPPER],
                                 text=True, capture_output=True, timeout=60)
        (args.output / "invalid.log").write_text(invalid.stdout + invalid.stderr)
        shutil.copytree(app / "out", args.output / "invalid-out", dirs_exist_ok=True)
        assert invalid.returncode != 0 and not (app / "out/prediction.json").exists()
        check_platform_outputs(app / "out", canaries, exit_code=invalid.returncode)
        print("PASS invalid structure: nonzero exit, no prediction JSON", flush=True)


if __name__ == "__main__":
    main()

