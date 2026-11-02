"""First-model structures and headless PyMOL sessions for ZincSight results.

Predicted Zn atoms are visualization markers, not assignments of metal identity.
The original atoms remain present, including experimentally observed ions.
"""
from __future__ import annotations

import copy
from io import StringIO
import json
import math
from pathlib import Path
import re
import string

CITATION = (
    "Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion "
    "location in proteins with ZincSight. Protein Science. 2025;34(11):e70350. "
    "https://doi.org/10.1002/pro.70350"
)
MARKER_NOTE = "Predicted Zn labels mark potential binding sites; they do not identify the bound metal."
CHAIN_CHARACTERS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def _validate_sites(result):
    sites = result.get("sites", [])
    if not isinstance(sites, list):
        raise ValueError("Structure export expects a list of predicted sites.")
    for site in sites:
        coordinates = site.get("coordinates_angstrom", [])
        probability = site.get("probability")
        if (len(coordinates) != 3 or not all(isinstance(value, (int, float))
                and not isinstance(value, bool) and math.isfinite(value) for value in coordinates)):
            raise ValueError("Every predicted marker must have three finite coordinates.")
        if (not isinstance(probability, (int, float)) or isinstance(probability, bool)
                or not math.isfinite(probability) or not 0 <= probability <= 1):
            raise ValueError("Every predicted marker must have a finite probability between zero and one.")
    return sites


def _first_model(structure_text, file_format):
    from Bio.PDB import MMCIFParser, PDBParser
    if file_format not in ("pdb", "cif", "mmcif"):
        raise ValueError("Structure export supports PDB and mmCIF input.")
    parser = PDBParser(QUIET=True) if file_format == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("input", StringIO(structure_text))
    model = next(structure.get_models(), None)
    if model is None or not list(model.get_atoms()):
        raise ValueError("The input has no atoms in its first structural model.")
    return model


def _first_model_cif(structure_text, file_format, model):
    """Keep original mmCIF coordinate strings rather than rounding via Atom objects."""
    from Bio.PDB import MMCIFIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    if file_format != "pdb":
        data = MMCIF2Dict(StringIO(structure_text))
        models = data.get("_atom_site.pdbx_PDB_model_num")
        if models:
            keep = [i for i, value in enumerate(models) if value == models[0]]
            for name, values in data.items():
                if name.startswith("_atom_site."):
                    data[name] = [values[i] for i in keep]
        return data
    writer = MMCIFIO()
    writer.set_structure(model)
    buffer = StringIO()
    writer.save(buffer)
    return MMCIF2Dict(StringIO(buffer.getvalue()))


def _write_cif(data, path):
    from Bio.PDB import MMCIFIO
    writer = MMCIFIO()
    writer.set_dict(data)
    writer.save(str(path))


def _marker_chain(data):
    used = set(data.get("_atom_site.auth_asym_id", [])) | set(data.get("_atom_site.label_asym_id", []))
    for candidate in "Z" + CHAIN_CHARACTERS:
        if candidate not in used:
            return candidate
    candidate = "PREDICTED"
    while candidate in used:
        candidate += "_ZN"
    return candidate


def _combined_cif(data, sites, marker_chain):
    combined = copy.deepcopy(data)
    columns = [name for name in combined if name.startswith("_atom_site.")]
    row_count = len(combined["_atom_site.id"])
    for name, default in {"_atom_site.auth_asym_id": "_atom_site.label_asym_id",
                          "_atom_site.auth_seq_id": "_atom_site.label_seq_id",
                          "_atom_site.auth_comp_id": "_atom_site.label_comp_id",
                          "_atom_site.auth_atom_id": "_atom_site.label_atom_id"}.items():
        if name not in combined:
            combined[name] = list(combined[default])
            columns.append(name)
    numeric_ids = [int(value) for value in combined["_atom_site.id"] if str(value).isdigit()]
    first_id = max(numeric_ids, default=row_count) + 1
    model_number = combined.get("_atom_site.pdbx_PDB_model_num", ["1"])[0]
    for rank, site in enumerate(sites, 1):
        x, y, z = site["coordinates_angstrom"]
        values = {
            "group_PDB": "HETATM", "id": str(first_id + rank - 1), "type_symbol": "ZN",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": marker_chain, "label_entity_id": ".", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": repr(x), "Cartn_y": repr(y), "Cartn_z": repr(z),
            "occupancy": "1.00", "B_iso_or_equiv": "0.00", "pdbx_formal_charge": "?",
            "auth_seq_id": str(rank), "auth_comp_id": "ZN", "auth_asym_id": marker_chain,
            "auth_atom_id": "ZN", "pdbx_PDB_model_num": model_number,
        }
        for name in columns:
            combined[name].append(values.get(name.removeprefix("_atom_site."), "?"))
    combined["_zincsight_annotation.method"] = "ZincSight"
    combined["_zincsight_annotation.citation"] = CITATION
    combined["_zincsight_annotation.marker_note"] = MARKER_NOTE
    combined["_zincsight_annotation.predicted_marker_chain"] = marker_chain
    combined["_zincsight_annotation.model_policy"] = "First structural model only; all original atoms retained."
    if sites:
        combined["_zincsight_site.id"] = [str(i) for i in range(1, len(sites) + 1)]
        combined["_zincsight_site.probability"] = [repr(site["probability"]) for site in sites]
        combined["_zincsight_site.coordinating_residues"] = [
            ",".join(site.get("coordinating_residues", [])) or "." for site in sites]
    return combined


def _write_pdb(model, sites, marker_chain, path, warnings):
    from Bio.PDB import PDBIO
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    import numpy as np

    if len(sites) > 9999:
        warnings.append("PDB omitted: predicted marker residue numbering exceeds the PDB format; use mmCIF/PyMOL.")
        return None, {}, {}, None
    export_model = copy.deepcopy(model)
    # PDBIO writes every disordered residue/alternate atom, whereas get_atoms()
    # exposes only the selected conformer. Validate exactly the atoms it writes.
    original_residues = [residue for chain in export_model.get_chains()
                         for residue in chain.get_unpacked_list()]
    original_atoms = [atom for residue in original_residues
                      for atom in residue.get_unpacked_list()]
    if (any(len(residue.get_resname()) > 3 for residue in original_residues)
            or any(len(atom.get_fullname().strip()) > 4 for atom in original_atoms)):
        warnings.append("PDB omitted: component or atom names exceed fixed-width PDB fields; use mmCIF/PyMOL.")
        return None, {}, {}, None
    if len(original_atoms) + len(sites) > 99999:
        warnings.append("PDB omitted: atom serial numbers exceed the PDB format; use mmCIF/PyMOL.")
        return None, {}, {}, None
    coordinates = [float(value) for atom in original_atoms for value in atom.coord]
    coordinates.extend(float(value) for site in sites for value in site["coordinates_angstrom"])
    if any(not math.isfinite(value) or not -999.999 <= value <= 9999.999 for value in coordinates):
        warnings.append("PDB omitted: coordinates exceed the fixed-width PDB format; use mmCIF/PyMOL.")
        return None, {}, {}, None
    original_chains = list(export_model.get_chains())
    available = list(CHAIN_CHARACTERS)
    if len(original_chains) + bool(sites) > len(available) + any(chain.id == " " for chain in original_chains):
        warnings.append("PDB omitted: this structure needs more than 62 chain identifiers; use the mmCIF and PyMOL files.")
        return None, {}, {}, None
    chain_mapping = {}
    for chain in original_chains:
        if chain.id == " ":
            chain_mapping[chain.id] = chain.id
        elif len(chain.id) == 1 and chain.id in available:
            chain_mapping[chain.id] = chain.id
            available.remove(chain.id)
    for chain in original_chains:
        if chain.id not in chain_mapping:
            chain_mapping[chain.id] = available.pop(0)
    pdb_marker_chain = marker_chain if marker_chain in available else available[0] if available else None

    # Detach before renaming to avoid collisions between old/new chain IDs.
    for chain in original_chains:
        export_model.detach_child(chain.id)
    residue_mapping = {}
    for chain in original_chains:
        original_chain = chain.id
        chain.id = chain_mapping[original_chain]
        residues = list(chain.get_residues())
        renumber = any(not -999 <= residue.id[1] <= 9999 for residue in residues)
        if renumber:
            if len(residues) > 9999:
                warnings.append("PDB omitted: residue numbering exceeds PDB limits; use mmCIF or PyMOL.")
                return None, chain_mapping, {}, None
            for residue in residues:
                chain.detach_child(residue.id)
            for number, residue in enumerate(residues, 1):
                original = f"{original_chain}_{residue.id[1]}{residue.id[2].strip()}"
                residue.id = (residue.id[0], number, " ")
                chain.add(residue)
                residue_mapping[original] = {"chain": chain.id, "resi": str(number)}
        export_model.add(chain)
    if any(old != new for old, new in chain_mapping.items()):
        warnings.append("PDB chain identifiers were mapped to single characters; original identifiers remain in mmCIF/PyMOL and the annotations JSON.")
    if residue_mapping:
        warnings.append("PDB residues were renumbered to fit the format; see annotations JSON for the mapping.")
    if sites:
        predicted = Chain(pdb_marker_chain)
        for rank, site in enumerate(sites, 1):
            residue = Residue(("H_ZN", rank, " "), "ZN", " ")
            residue.add(Atom("ZN", np.asarray(site["coordinates_angstrom"], dtype=float),
                             0.0, 1.0, " ", " ZN ", rank, element="ZN"))
            predicted.add(residue)
        export_model.add(predicted)
    writer = PDBIO()
    writer.set_structure(export_model)
    buffer = StringIO()
    try:
        writer.save(buffer)
    except (ValueError, TypeError) as exc:
        warnings.append(f"PDB omitted because the input exceeds PDB formatting limits: {exc}. Use mmCIF/PyMOL.")
        return None, chain_mapping, residue_mapping, None
    remarks = [
        "REMARK 900 ZINCSIGHT PREDICTED METAL-BINDING SITE MARKERS",
        "REMARK 900 DOI: 10.1002/pro.70350",
        "REMARK 900 ZN MARKERS ARE PREDICTIONS, NOT ASSIGNMENTS OF METAL IDENTITY.",
        "REMARK 900 FIRST MODEL ONLY; ALL ORIGINAL ATOMS AND EXPERIMENTAL IONS RETAINED.",
        f"REMARK 900 PREDICTED MARKER CHAIN: {pdb_marker_chain or 'NONE'}; RESIDUE = SITE RANK.",
        "REMARK 900 MARKER OCCUPANCY=1 AND B=0 ARE DISPLAY PLACEHOLDERS, NOT MEASUREMENTS.",
    ]
    for old, new in chain_mapping.items():
        if old != new:
            remarks.append(f"REMARK 900 PDB CHAIN MAP {old!r} -> {new!r}; SEE ANNOTATIONS JSON.")
    for rank, site in enumerate(sites, 1):
        remarks.append(f"REMARK 910 SITE {rank} PROBABILITY {site['probability']:.8g}")
        remarks.append(f"REMARK 910 SITE {rank} RESIDUES {','.join(site.get('coordinating_residues', []))}")
    path.write_text("\n".join(line[:80] for line in remarks) + "\n" + buffer.getvalue(), encoding="utf-8")
    return path, chain_mapping, residue_mapping, pdb_marker_chain


def _pymol_session(original_cif, sites, marker_chain, path, warnings):
    import pymol2
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.load(str(original_cif), "original_structure", state=1)
        cmd.hide("everything", "all")
        cmd.bg_color("white")
        cmd.set_color("zincsight_protein", [0.85, 0.85, 1.0])
        cmd.set_color("zincsight_binding_carbon", [0.6, 0.4, 0.8])
        cmd.color("zincsight_protein", "original_structure and polymer.protein")
        cmd.show("cartoon", "original_structure and polymer.protein")
        cmd.show("cartoon", "original_structure and polymer.nucleic")
        cmd.color("orange", "original_structure and polymer.nucleic")
        cmd.show("sticks", "original_structure and not polymer and not solvent")
        cmd.select("experimental_metals", "original_structure and elem Zn+Fe+Mn+Co+Cu+Ni+Mg+Ca")
        cmd.show("spheres", "experimental_metals")
        cmd.set("sphere_scale", 0.45, "experimental_metals")
        original_atoms = cmd.get_model("original_structure").atom
        all_binding = set()
        for rank, site in enumerate(sites, 1):
            selected = set()
            for identifier in site.get("coordinating_residues", []):
                chain, separator, resi = str(identifier).rpartition("_")
                if not separator:
                    warnings.append(f"Site {rank}: could not interpret coordinating residue {identifier!r}.")
                    continue
                matched = [atom for atom in original_atoms if
                           (atom.chain == chain or (not atom.chain.strip() and not chain.strip())) and
                           (atom.resi == resi or (re.fullmatch(r"-?\d+", resi) and
                            re.match(r"-?\d+", atom.resi) and re.match(r"-?\d+", atom.resi).group() == resi))]
                if not matched:
                    warnings.append(f"Site {rank}: coordinating residue {identifier!r} was not found for highlighting.")
                if len({atom.resi for atom in matched}) > 1:
                    warnings.append(f"Site {rank}: insertion-code ambiguity for {identifier!r}; all matching residues highlighted.")
                selected.update(atom.index for atom in matched)
            selection = f"site_{rank:03d}_coordinating"
            expression = "original_structure and index " + "+".join(map(str, sorted(selected))) if selected else "none"
            cmd.select(selection, expression)
            cmd.disable(selection)
            all_binding.update(selected)
            cmd.pseudoatom("predicted_sites", pos=list(site["coordinates_angstrom"]), elem="Zn",
                           name="ZN", resn="ZN", chain=marker_chain, resi=str(rank), vdw=1.1)
            marker = f"predicted_sites and resi {rank}"
            probability = float(site["probability"])
            color = f"zincsight_probability_{rank:03d}"
            cmd.set_color(color, [1.0 - probability, 0.15 + 0.65 * probability, 0.15])
            cmd.color(color, marker)
            cmd.label(marker, repr(f"Predicted site {rank} | p={probability:.3g}"))
        if all_binding:
            expression = "original_structure and index " + "+".join(map(str, sorted(all_binding)))
            cmd.select("coordinating_residues", expression)
            cmd.show("sticks", "coordinating_residues")
            cmd.color("zincsight_binding_carbon", "coordinating_residues and elem C")
            for element, color in [("N", "blue"), ("O", "red"), ("S", "yellow")]:
                cmd.color(color, f"coordinating_residues and elem {element}")
            cmd.disable("coordinating_residues")
        if sites:
            cmd.show("spheres", "predicted_sites")
            cmd.show("labels", "predicted_sites")
            cmd.set("sphere_scale", 0.65, "predicted_sites")
            cmd.set("label_color", "black", "predicted_sites")
            cmd.set("label_outline_color", "white", "predicted_sites")
            cmd.set("label_size", 14, "predicted_sites")
        cmd.disable("experimental_metals")
        cmd.set_title("original_structure", 1, "ZincSight | Protein Science (2025) | doi:10.1002/pro.70350")
        cmd.orient("original_structure")
        cmd.zoom("all", buffer=3)
        cmd.save(str(path))


def export_structure_outputs(structure_text, file_format, result, output_dir, stem="structure"):
    """Export first-model input atoms, predicted markers, annotations, and a session.

    PDB is rounded to its standard 0.001 Å coordinate precision. mmCIF retains
    the input's original coordinate strings and the full predicted coordinates.
    Unsupported PDB chain/residue limits produce a warning and CIF/PSE fallback.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,99}", stem) or ".." in stem:
        raise ValueError("Use a safe structure output stem without directory components.")
    sites = _validate_sites(result)
    model = _first_model(structure_text, file_format)
    original_data = _first_model_cif(structure_text, file_format, model)
    marker_chain = _marker_chain(original_data)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    cif_path, pdb_path = output / f"{stem}.cif", output / f"{stem}.pdb"
    pse_path, annotations_path = output / f"{stem}.pse", output / f"{stem}_sites.json"
    warnings = []
    _write_cif(_combined_cif(original_data, sites, marker_chain), cif_path)
    pdb_path, chain_mapping, residue_mapping, pdb_marker_chain = _write_pdb(
        model, sites, marker_chain, pdb_path, warnings)
    # A temporary original-only CIF lets the session distinguish observed atoms
    # from computational markers without deleting or reclassifying input ions.
    import tempfile
    with tempfile.TemporaryDirectory(prefix="zincsight-export-") as temporary:
        original_cif = Path(temporary) / "original.cif"
        _write_cif(original_data, original_cif)
        _pymol_session(original_cif, sites, marker_chain, pse_path, warnings)
    annotations = {
        "method": "ZincSight", "citation": CITATION, "marker_note": MARKER_NOTE,
        "model_policy": "First structural model only; all original atoms retained.",
        "coordinate_frame": "Input structure coordinate frame; no alignment or recentering.",
        "pdb_coordinate_precision_angstrom": 0.001,
        "predicted_marker_chain": marker_chain, "pdb_predicted_marker_chain": pdb_marker_chain,
        "predicted_marker_residue_name": "ZN", "predicted_marker_element": "Zn",
        "pdb_chain_mapping": chain_mapping, "pdb_residue_mapping": residue_mapping,
        "visualization": {"original_structure": "lavender protein cartoon; original ions retained separately",
                          "coordinating_residues": "purple carbon sticks with element-colored donor atoms",
                          "predicted_sites": "Zn-element spheres, colored from red (low probability) to green (high probability), with probability labels"},
        "site_count": len(sites), "sites": copy.deepcopy(sites), "warnings": warnings,
    }
    annotations_path.write_text(json.dumps(annotations, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    artifacts = {"mmcif": str(cif_path), "pymol_session": str(pse_path), "site_annotations": str(annotations_path)}
    if pdb_path:
        artifacts["pdb"] = str(pdb_path)
    return {"artifact_paths": artifacts, **{key: annotations[key] for key in (
        "pdb_chain_mapping", "pdb_residue_mapping", "predicted_marker_chain",
        "pdb_predicted_marker_chain", "predicted_marker_residue_name", "warnings")}}
