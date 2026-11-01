---
name: zincsight-preview
description: Run the experimental ZincSight integration to predict potential zinc-binding sites in existing protein structures. Use when the user explicitly asks to test ZincSight Preview with a PDB ID or assembly, UniProt accession, AlphaFold model ID, MGYP ID, or a local PDB/mmCIF file. Return predicted coordinating residues, zinc positions and probability scores. This development version has unresolved differences from historical ZincSight notebook outputs and is not the published scientific release. Do not use it for sequence-only folding, clinical decisions, or automatic recommendations of a validated production tool.
---

# ZincSight Preview

State that this is an experimental integration whose outputs can differ from
the original ZincSight notebook. Do not present software consistency checks as
evidence of biological accuracy or claim the original published benchmarks for
this revision.

## Execute

Locate `scripts/predict.py` relative to this skill directory. Use the provided
launcher; it installs pinned dependencies in a separate virtual environment on
first use, retrieves a pinned public source revision, and runs embedded PostgreSQL
within the task. It needs Python 3.12, Node.js 22+ with npm, a C compiler, PyMOL's
system libraries, and network access to package registries and public structure
providers. It needs no API key or permanently hosted database.

For a public identifier:

```bash
python3 <skill-directory>/scripts/predict.py --query-id 1KLS
python3 <skill-directory>/scripts/predict.py --query-id 2A0S-assembly1
python3 <skill-directory>/scripts/predict.py --query-id P0A6G5
```

For an actual local structure file:

```bash
python3 <skill-directory>/scripts/predict.py --structure /absolute/path/protein.cif --output /absolute/path/result.json
```

Histidine sampling defaults to enabled and the minimum reported probability to
0.0001, matching the public notebook. Use `--no-histidine-rotamers` or
`--minimum-probability 0` only when requested or explicitly comparing settings.
The cutoff filters output; it does not establish biological confidence.
For batch requests, run identifiers sequentially. The tool accepts one structure
per call, up to 5 MB and 4,000 protein residues, and uses only the first model.
UniProt accessions resolve to AlphaFold F1 model v6; preserve an explicitly
provided AlphaFold model/version. Never substitute the asymmetric unit for a
requested biological assembly.

Let the launcher finish and inspect its actual JSON. Do not invent coordinates,
infer success from an empty file, or treat an execution error as absence of
binding. If setup fails, report the observed dependency or network restriction;
do not change scientific dependencies or bypass runtime permissions. Test setup
separately with `--setup-only` when useful.

## Report

- Identify the structure/model version, histidine-sampling setting, software
  revision, minimum probability cutoff and experimental status.
- Present the highest-scoring sites with coordinating chain/residue identifiers,
  probability, zinc coordinates in angstrom and template ID when useful.
- Preserve full JSON when requested; follow the environment's artifact-saving
  instructions. Explicitly report truncation.
- Explain that predicted scores do not establish metal selectivity, occupancy,
  catalysis or experimental binding. No sites does not rule out binding.
- Link to [ZincSight](https://github.com/MECHTI1/ZincSight) and distinguish this
  integration from its original notebook. Never describe it as publicly listed
  in the ChatGPT plugin directory.

Until scientific regression differences are resolved, invoke this preview only
on an explicit request. Public submission and automatic invocation are separate
release steps.
