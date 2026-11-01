"""Executed only inside a temporary copy of the predictor."""
import json
import os
from pathlib import Path
import sys


def write_first_model(model, input_path, file_format, output_path):
    """Keep mmCIF coordinate text intact while applying the first-model policy."""
    from Bio.PDB import MMCIFIO
    writer = MMCIFIO()
    if file_format == 'cif':
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        data = MMCIF2Dict(str(input_path))
        model_numbers = data.get('_atom_site.pdbx_PDB_model_num')
        if model_numbers and len(set(model_numbers)) > 1:
            selected = [i for i, value in enumerate(model_numbers) if value == model_numbers[0]]
            for key, values in data.items():
                if key.startswith('_atom_site.'):
                    data[key] = [values[i] for i in selected]
        # Serializing Bio.PDB Atom objects rounds coordinates to 3 decimals.
        # Re-emit original mmCIF fields instead, preserving higher precision.
        writer.set_dict(data)
    else:
        writer.set_structure(model)
    writer.save(str(output_path))


def main():
    from Bio.PDB import MMCIFParser, PDBParser, MMCIFIO

    root = Path(__file__).resolve().parent
    file_format, rotamers = sys.argv[1:]
    parser = PDBParser(QUIET=True) if file_format == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("query", str(root / f"input.{file_format}"))
    model = next(structure.get_models())
    residues = [r for r in model.get_residues() if r.id[0] == " "]
    max_residues = int(os.environ.get('ZINCSIGHT_MAX_RESIDUES', '4000'))
    if not 3 <= len(residues) <= max_residues:
        raise ValueError(f"Input must have between 3 and {max_residues} protein residues in its first model.")
    query = root / "query"
    query.mkdir()
    write_first_model(model, root / f'input.{file_format}', file_format, query / 'query.cif')
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
