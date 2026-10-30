"""Executed only inside a temporary copy of the predictor."""
import json
import os
from pathlib import Path
import sys


def main():
    from Bio.PDB import MMCIFParser, PDBParser, MMCIFIO

    root = Path(__file__).resolve().parent
    file_format, rotamers = sys.argv[1:]
    parser = PDBParser(QUIET=True) if file_format == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("query", str(root / f"input.{file_format}"))
    model = next(structure.get_models())
    residues = [r for r in model.get_residues() if r.id[0] == " "]
    if not 3 <= len(residues) <= 2000:
        raise ValueError("Input must have between 3 and 2000 protein residues in its first model.")
    query = root / "query"
    query.mkdir()
    writer = MMCIFIO()
    writer.set_structure(model)
    writer.save(str(query / "query.cif"))
    if os.environ.get("ZINCSIGHT_DATABASE_BACKEND", "native") == "pglite":
        from embedded_database import embedded_database
        with embedded_database(root, Path(os.environ["ZINCSIGHT_PGLITE_LAUNCHER"])):
            from src.primary_script import main as predict
            # Input is already validated and converted to mmCIF. The bundled
            # database has been loaded by the adapter, so enter the same pipeline
            # after the standard entry point's download/psql preparation steps.
            result = predict([str(query / "query.cif")], rotamers == "1", str(root / "results"), 2, False)
    else:
        from main_execute import execute_zincsight
        result = execute_zincsight(rotamers == "1", "", str(query), str(root / "results"), 2, False)
    if result is not False and not (root / "results/table/table_predicted_zn_binding_sites.csv").is_file():
        raise RuntimeError("Missing prediction output.")
    (root / "status.json").write_text(json.dumps({"has_sites": result is not False}))


if __name__ == "__main__":
    main()
