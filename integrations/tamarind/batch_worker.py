"""Disposable batch process; scientific source and normalization stay pinned."""
import json
import os
from pathlib import Path
import sys


def normalize_structure(item, root, query, max_residues, write_first_model):
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import is_aa

    try:
        parser = PDBParser(QUIET=True) if item["file_format"] == "pdb" else MMCIFParser(QUIET=True)
        original = root / item["path"]
        structure = parser.get_structure(item["key"], str(original))
        model = next(structure.get_models())
        count = sum(residue.id[0] == " " and is_aa(residue, standard=True)
                    for residue in model.get_residues())
        if not 3 <= count <= max_residues:
            raise ValueError(f"Input must have between 3 and {max_residues} standard protein residues in its first model; found {count}.")
        normalized = query / f"{item['key']}.cif"
        write_first_model(model, original, item["file_format"], normalized)
        return normalized, count
    except Exception as exc:
        raise ValueError(f"Structure {item['key']} ({item['source']}): {exc}") from exc


def main():
    from embedded_database import embedded_database
    from pinned_worker import write_first_model

    root = Path(__file__).resolve().parent
    rotamers, cpu_workers = sys.argv[1:]
    manifest = json.loads((root / "batch-inputs.json").read_text(encoding="utf-8"))
    query = root / "query"
    query.mkdir()
    maximum = int(os.environ.get("ZINCSIGHT_MAX_RESIDUES", "4000"))
    normalized, counts = [], {}
    for item in manifest:
        path, count = normalize_structure(item, root, query, maximum, write_first_model)
        normalized.append(str(path))
        counts[item["key"]] = count
    # All originals are validated before starting one private database. Import
    # after connection installation: scientific modules import it by name.
    with embedded_database(root, Path(os.environ["ZINCSIGHT_PGLITE_LAUNCHER"])):
        from src.primary_script import main as predict
        result = predict(normalized, rotamers == "1", str(root / "results"), int(cpu_workers), False)
    if result is not False and not (root / "results/table/table_predicted_zn_binding_sites.csv").is_file():
        raise RuntimeError("Missing batch prediction output.")
    (root / "status.json").write_text(json.dumps({"has_sites": result is not False,
                                                "protein_residue_counts": counts}), encoding="utf-8")


if __name__ == "__main__":
    main()
