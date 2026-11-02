"""Structural and PyMOL-session roundtrip checks for exported predictions."""
import copy
from io import StringIO
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pymol2

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "structure_outputs", ROOT / "integrations/tamarind/structure_outputs.py")
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)

PDB = """MODEL        1
ATOM      1  N   CYS A   1       1.000   2.000   3.000  1.00 20.00           N
ATOM      2  CA  CYS A   1       2.000   3.000   4.000  1.00 20.00           C
ATOM      3  SG  CYS A   1       3.000   4.000   5.000  1.00 20.00           S
ATOM      4  ND1 HIS A   2       4.000   5.000   6.000  1.00 20.00           N
HETATM    5 ZN    ZN Z   1       7.000   8.000   9.000  1.00 10.00          ZN
ENDMDL
MODEL        2
ATOM      6  CA  CYS A   1      92.000  93.000  94.000  1.00 20.00           C
ENDMDL
END
"""
SITE = {
    "probability": 0.91234567, "coordinates_angstrom": [3.123456789, 4.23456789, 5.3456789],
    "coordinating_residues": ["A_1", "A_2"], "residue_combination": "{C,H}",
    "template_id": "19", "template_rmsd_angstrom": 0.25,
}


class StructureOutputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.result = {"site_count": 1, "sites": [copy.deepcopy(SITE)]}

    def export(self, content=PDB, file_format="pdb"):
        return exporter.export_structure_outputs(content, file_format, self.result, self.output)

    def test_first_model_original_ion_and_marker_survive_file_and_session_roundtrip(self):
        before = copy.deepcopy(self.result)
        exported = self.export()
        paths = exported["artifact_paths"]
        parsed = PDBParser(QUIET=True).get_structure("export", paths["pdb"])
        self.assertEqual(len(list(parsed.get_models())), 1)
        model = next(parsed.get_models())
        self.assertEqual(len(list(model.get_atoms())), 6)
        self.assertEqual(tuple(next(model["Z"].get_atoms()).coord), (7.0, 8.0, 9.0))
        predicted = next(model[exported["pdb_predicted_marker_chain"]].get_atoms())
        for actual, expected in zip(predicted.coord, SITE["coordinates_angstrom"]):
            self.assertAlmostEqual(float(actual), expected, places=3)
        text = Path(paths["pdb"]).read_text()
        self.assertIn("10.1002/pro.70350", text)
        self.assertIn("NOT ASSIGNMENTS OF METAL IDENTITY", text)
        annotations = json.loads(Path(paths["site_annotations"]).read_text())
        self.assertEqual(annotations["sites"], [SITE])
        self.assertEqual(before, self.result)
        with pymol2.PyMOL() as session:
            cmd = session.cmd
            cmd.load(paths["pymol_session"])
            self.assertEqual(cmd.count_states("original_structure"), 1)
            self.assertEqual(cmd.count_atoms("original_structure"), 5)
            self.assertEqual(cmd.count_atoms("experimental_metals"), 1)
            self.assertEqual(cmd.count_atoms("predicted_sites"), 1)
            self.assertEqual(cmd.count_atoms("coordinating_residues"), 4)
            atom = cmd.get_model("predicted_sites").atom[0]
            for actual, expected in zip(atom.coord, SITE["coordinates_angstrom"]):
                self.assertAlmostEqual(actual, expected, places=5)
            labels = []
            cmd.iterate("predicted_sites", "labels.append(label)", space={"labels": labels})
            self.assertIn("Predicted site 1", labels[0])
            self.assertIn("p=0.912", labels[0])

    def test_mmcif_precision_original_chain_ids_and_pdb_mapping_are_explicit(self):
        model = next(PDBParser(QUIET=True).get_structure("input", StringIO(PDB)).get_models())
        data = exporter._first_model_cif(PDB, "pdb", model)
        data["_atom_site.Cartn_x"][0] = "1.123456789012345"
        for i, chain in enumerate(data["_atom_site.auth_asym_id"]):
            if chain == "A":
                data["_atom_site.auth_asym_id"][i] = "LONG_CHAIN"
                if data["_atom_site.auth_seq_id"][i] == "1":
                    data["_atom_site.auth_seq_id"][i] = "10001"
        writer = MMCIFIO()
        writer.set_dict(data)
        buffer = StringIO()
        writer.save(buffer)
        self.result["sites"][0]["coordinating_residues"] = ["LONG_CHAIN_10001", "LONG_CHAIN_2"]
        exported = self.export(buffer.getvalue(), "cif")
        cif = MMCIF2Dict(exported["artifact_paths"]["mmcif"])
        self.assertEqual(cif["_atom_site.Cartn_x"][0], "1.123456789012345")
        self.assertIn("LONG_CHAIN", cif["_atom_site.auth_asym_id"])
        self.assertEqual(cif["_atom_site.Cartn_x"][-1], repr(SITE["coordinates_angstrom"][0]))
        self.assertNotEqual(exported["pdb_chain_mapping"]["LONG_CHAIN"], "LONG_CHAIN")
        self.assertIn("LONG_CHAIN_10001", exported["pdb_residue_mapping"])
        self.assertTrue(exported["warnings"])
        with pymol2.PyMOL() as session:
            session.cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertIn("LONG_CHAIN", session.cmd.get_chains("original_structure"))
            self.assertEqual(session.cmd.count_atoms("coordinating_residues"), 4)

    def test_empty_prediction_keeps_structure_without_inventing_metal_markers(self):
        self.result = {"site_count": 0, "sites": []}
        exported = self.export()
        parsed = PDBParser(QUIET=True).get_structure("export", exported["artifact_paths"]["pdb"])
        self.assertEqual(len(list(parsed.get_atoms())), 5)
        with pymol2.PyMOL() as session:
            session.cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertNotIn("predicted_sites", session.cmd.get_names("objects"))
            self.assertEqual(session.cmd.count_atoms("experimental_metals"), 1)
        annotations = json.loads(Path(exported["artifact_paths"]["site_annotations"]).read_text())
        self.assertEqual(annotations["site_count"], 0)

    def test_blank_chain_identifier_is_preserved_and_highlighted(self):
        lines = []
        for line in PDB.splitlines():
            if line.startswith("ATOM  ") and line[21] == "A":
                line = line[:21] + " " + line[22:]
            lines.append(line)
        self.result["sites"][0]["coordinating_residues"] = ["_1", "_2"]
        exported = self.export("\n".join(lines) + "\n")
        self.assertEqual(exported["pdb_chain_mapping"][" "], " ")
        with pymol2.PyMOL() as session:
            session.cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertEqual(session.cmd.count_atoms("coordinating_residues"), 4)

    def _cif_data(self):
        model = next(PDBParser(QUIET=True).get_structure("input", StringIO(PDB)).get_models())
        return exporter._first_model_cif(PDB, "pdb", model)

    def _export_cif_fallback(self, data):
        writer = MMCIFIO()
        writer.set_dict(data)
        text = StringIO()
        writer.save(text)
        exported = self.export(text.getvalue(), "cif")
        self.assertNotIn("pdb", exported["artifact_paths"])
        self.assertFalse((self.output / "structure.pdb").exists())
        with pymol2.PyMOL() as session:
            session.cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertEqual(session.cmd.count_atoms("predicted_sites"), 1)
        return exported, MMCIF2Dict(exported["artifact_paths"]["mmcif"])

    def test_long_component_name_uses_cif_and_pymol_instead_of_shifted_pdb_fields(self):
        data = self._cif_data()
        data["_atom_site.label_comp_id"][-1] = "LONG"
        exported, saved = self._export_cif_fallback(data)
        self.assertEqual(saved["_atom_site.label_comp_id"][4], "LONG")
        self.assertTrue(any("names exceed" in warning for warning in exported["warnings"]))

    def test_long_atom_name_uses_cif_and_pymol_instead_of_shifted_pdb_fields(self):
        data = self._cif_data()
        data["_atom_site.label_atom_id"][0] = "LONGATOM"
        exported, saved = self._export_cif_fallback(data)
        self.assertEqual(saved["_atom_site.label_atom_id"][0], "LONGATOM")
        self.assertTrue(any("names exceed" in warning for warning in exported["warnings"]))

    def test_unselected_altloc_coordinate_overflow_uses_cif_and_pymol(self):
        data = self._cif_data()
        # The selected A conformer fits PDB; PDBIO also emits B, whose X value
        # cannot fit its 8.3 fixed-width field. Neither original may be dropped.
        data["_atom_site.label_alt_id"][0] = "A"
        data["_atom_site.occupancy"][0] = "0.75"
        for name, values in data.items():
            if name.startswith("_atom_site."):
                values.append(values[0])
        data["_atom_site.id"][-1] = "99"
        data["_atom_site.label_alt_id"][-1] = "B"
        data["_atom_site.occupancy"][-1] = "0.25"
        data["_atom_site.Cartn_x"][-1] = "10000.000"
        exported, saved = self._export_cif_fallback(data)
        self.assertIn("10000.000", saved["_atom_site.Cartn_x"])
        self.assertEqual(len(saved["_atom_site.id"]), 7)
        self.assertTrue(any("coordinates exceed" in warning for warning in exported["warnings"]))
        with pymol2.PyMOL() as session:
            session.cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertEqual(session.cmd.count_atoms("original_structure"), 6)

    def test_unsafe_output_paths_and_nonfinite_markers_are_rejected(self):
        with self.assertRaises(ValueError):
            exporter.export_structure_outputs(PDB, "pdb", self.result, self.output, "../escape")
        self.result["sites"][0]["coordinates_angstrom"][0] = float("nan")
        with self.assertRaises(ValueError):
            self.export()
        self.assertEqual(list(self.output.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
