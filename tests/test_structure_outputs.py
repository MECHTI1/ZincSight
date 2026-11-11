"""Structural and PyMOL-session roundtrip checks for exported predictions."""
import copy
from io import StringIO
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

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

    def test_pymol_scenes_reduce_label_overlap_without_hiding_or_changing_sites(self):
        # Deliberately put the highest probability last: the overview must use
        # the score, not assume the first result is always the best one.
        self.result = {"site_count": 3, "sites": [
            {**copy.deepcopy(SITE), "probability": probability,
             "coordinates_angstrom": [3.1 + rank, 4.2, 5.3]}
            for rank, probability in enumerate((0.4, 0.7, 0.95))
        ]}
        before = copy.deepcopy(self.result)
        exported = self.export()
        with pymol2.PyMOL() as session:
            cmd = session.cmd
            cmd.load(exported["artifact_paths"]["pymol_session"])
            self.assertEqual(cmd.get_scene_list(), ["Overview", "All_site_labels"])
            self.assertEqual(cmd.get_setting_int("float_labels"), 1)
            self.assertEqual(cmd.get_setting_float("ray_label_specular"), 0)
            self.assertEqual(cmd.count_atoms("predicted_sites", state=1), 3)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep spheres", state=1), 3)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep labels", state=1), 1)
            self.assertEqual(cmd.get_model("predicted_sites and rep labels", state=1).atom[0].resi, "3")
            label_positions = []
            cmd.iterate("predicted_sites and rep labels", "label_positions.append(s.label_position)",
                        space={"label_positions": label_positions})
            self.assertEqual(label_positions, [(0.0, 15.0, 20.0)])
            originals = cmd.get_model("original_structure", state=1).atom
            markers = cmd.get_model("predicted_sites", state=1).atom
            original_snapshot = [(atom.chain, atom.resi, atom.name, atom.symbol, tuple(atom.coord))
                                 for atom in originals]
            marker_snapshot = [(atom.chain, atom.resi, atom.name, atom.symbol, tuple(atom.coord))
                               for atom in markers]
            camera = cmd.get_view()
            labels = []
            cmd.iterate("predicted_sites", "labels.append(label)", space={"labels": labels})
            self.assertEqual(len(labels), 3)
            for rank, (label, site, atom) in enumerate(zip(labels, self.result["sites"], markers), 1):
                self.assertEqual(label, f"Predicted site {rank} | p={site['probability']:.3g}")
                for actual, expected in zip(atom.coord, site["coordinates_angstrom"]):
                    self.assertAlmostEqual(actual, expected, places=5)
            cmd.scene("All_site_labels", "recall", animate=0)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep labels", state=1), 3)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep spheres", state=1), 3)
            self.assertEqual(cmd.get_view(), camera)
            cmd.scene("Overview", "recall", animate=0)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep labels", state=1), 1)
            self.assertEqual(cmd.count_atoms("predicted_sites and rep spheres", state=1), 3)
            self.assertEqual(cmd.get_view(), camera)
            self.assertEqual([(a.chain, a.resi, a.name, a.symbol, tuple(a.coord)) for a in
                              cmd.get_model("original_structure", state=1).atom], original_snapshot)
            self.assertEqual([(a.chain, a.resi, a.name, a.symbol, tuple(a.coord)) for a in
                              cmd.get_model("predicted_sites", state=1).atom], marker_snapshot)
        self.assertEqual(self.result, before)

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
            self.assertEqual(session.cmd.get_scene_list(), ["Overview", "All_site_labels"])
            for scene in ("All_site_labels", "Overview"):
                session.cmd.scene(scene, "recall", animate=0)
                self.assertNotIn("predicted_sites", session.cmd.get_names("objects"))
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

    def test_more_than_9999_sites_all_survive_in_cif_annotations_and_session_handoff(self):
        sites = [{**SITE, "probability": (index + 1) / 10001,
                  "coordinates_angstrom": [index / 10, (index % 91) / 11, -(index % 101) / 13],
                  "template_id": str(index + 1)} for index in range(10000)]
        self.result = {"site_count": len(sites), "sites": sites}
        def write_session(original_cif, received_sites, marker_chain, path, warnings):
            self.assertTrue(original_cif.is_file())
            self.assertIs(received_sites, sites)
            self.assertEqual(len(received_sites), 10000)
            path.write_bytes(b"mock session retains all supplied site records")
        # Session styling is exercised with real PyMOL above. This test checks
        # the large-result format boundary without generating 10,000 labels.
        with patch.object(exporter, "_pymol_session", side_effect=write_session) as session:
            exported = self.export()
        session.assert_called_once()
        self.assertNotIn("pdb", exported["artifact_paths"])
        self.assertFalse((self.output / "structure.pdb").exists())
        self.assertTrue(any("marker residue numbering" in warning for warning in exported["warnings"]))
        cif = MMCIF2Dict(exported["artifact_paths"]["mmcif"])
        marker_indexes = [i for i, chain in enumerate(cif["_atom_site.auth_asym_id"])
                          if chain == exported["predicted_marker_chain"]]
        self.assertEqual(len(marker_indexes), 10000)
        self.assertEqual(len(cif["_atom_site.id"]), 10005)
        self.assertEqual(cif["_atom_site.auth_seq_id"][marker_indexes[-1]], "10000")
        self.assertEqual(cif["_zincsight_site.id"], [str(i) for i in range(1, 10001)])
        for rank in (0, 9998, 9999):
            index = marker_indexes[rank]
            self.assertEqual([float(cif["_atom_site.Cartn_" + axis][index]) for axis in "xyz"],
                             sites[rank]["coordinates_angstrom"])
        annotations = json.loads(Path(exported["artifact_paths"]["site_annotations"]).read_text())
        self.assertEqual(annotations["site_count"], 10000)
        self.assertEqual(annotations["sites"], sites)
        self.assertTrue(Path(exported["artifact_paths"]["pymol_session"]).is_file())


if __name__ == "__main__":
    unittest.main()
