# ZincSight on Tamarind

This branch packages **ZincSight Preview v0.2.1** as a Tamarind Custom Tool.
Upload multiple protein structures or a structure archive in one job. Each input
gets a viewable structure, prediction tables, and a downloadable PyMOL session,
including inputs with no predicted sites.

Please cite **Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction
of zinc ion location in proteins with ZincSight. Protein Science.
2025;34(11):e70350.** [https://doi.org/10.1002/pro.70350](https://doi.org/10.1002/pro.70350).
The tool's Paper link, result JSON, structure annotations, and downloadable
CITATION.txt and the bundled CITATION.bib point to this paper. Tamarind manages its native
Citations tab separately; the current Custom Tool schema has no citation override,
so that tab may continue to show Tamarind's default reference.

## Connect and deploy

| Tamarind setting | Value |
| --- | --- |
| Repository | MECHTI1/ZincSight |
| Branch | feature/tamarind-integration |
| Source directory | Repository root |
| Suggested tool ID | zincsight-preview |
| Display name | ZincSight Preview |

1. Open [Tamarind Custom Tools](https://app.tamarind.bio/app/custom-tools).
   Link this repository and branch to the organization's existing tool, or create
   a tool if needed. Tamarind's GitHub App must have repository access.
2. Build/deploy the source with automatic publishing off during validation.
   Dependencies and the pinned engine are installed under /opt because Tamarind
   mounts the submitted repository over /app at runtime.
3. Run a test with a small public structure such as
   [1KLS.cif](https://files.rcsb.org/download/1KLS.cif). Check the Summary
   preview, the downloaded structure and .pse files, and the citation files.
   Also test a batch with a structure that has no predicted sites.
4. Publish the validated version to the organization and choose its visibility.
   A Tamarind test or production job consumes account compute.

See [Tamarind's Custom Tool documentation](https://docs.tamarind.bio/tamarind/custom-tools).
Custom tools are organization tools; public catalog inclusion is a separate
arrangement with Tamarind. The master, sqlite, and MCP branches do not require
edits to use this integration branch.

## Uploads and ID preparation

The tool page accepts either or both of these inputs:

| Input | Accepted content |
| --- | --- |
| structures | One or more uploaded .pdb, .cif, or .mmcif files |
| batch_archive | A ZIP, TAR, TAR.GZ, or TGZ containing PDB/mmCIF/ENT files, optionally gzip compressed; or a prepared ID ZIP |

The combined batch limit is **100 structures**, **5 MB per uncompressed
structure**, and **100 MB total uncompressed structure data**. Each first model
must contain 3–4,000 standard protein residues with the side-chain coordinates
needed for prediction. All protein chains in that model are analyzed. Later
models are not analyzed or exported. Archives are read without extracting their
paths; unsafe paths, links, duplicate archive member names, and oversized content
are rejected. Separate input files with the same basename retain distinct output
IDs.

Tamarind prediction jobs have no internet access. For a browser form with
comma/newline ID entry and optional file uploads, open
[Prepare a Tamarind batch in Colab](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/feature/tamarind-integration/integrations/tamarind/Prepare_Tamarind_Batch.ipynb).
Run its input form, prepare the ZIP, download it, and upload it in batch_archive.
The notebook embeds the preparer source and requires no API key or GPU;
predictions still run on Tamarind.

Alternatively, run the same bundled helper on a networked computer:

~~~bash
# Run from the repository root; Python's standard library is sufficient.
python3 integrations/tamarind/prepare_inputs.py \
  --ids "1KLS, 2A0S-assembly1, P12345" \
  --output zincsight_inputs.zip

# ID lists and local files/folders can be combined in one batch.
python3 integrations/tamarind/prepare_inputs.py \
  --ids-file ids.csv structures/ another_structure.pdb.gz \
  --output zincsight_inputs.zip
~~~

The helper accepts comma/whitespace ID text and .txt/.csv/.tsv lists. Supported
IDs are PDB entries and assemblies, UniProt accessions, complete AlphaFold model
IDs, MGYP identifiers, and complete TED domain IDs. UniProt accessions resolve
to AlphaFold F1 model v6; explicit AlphaFold fragments/versions and TED model
versions are preserved. Examples of explicit forms are AF-P12345-F2-model_v4 and
AF-A0A287EUI1-F1-model_v4_TED01. Availability depends on the public provider.

The helper also accepts local structures, folders, and ZIP/TAR archives. It
downloads only from an exact allowlist of public RCSB, AlphaFold, ESM Atlas,
and TED endpoints, with bounded download sizes and retries. It never reads
account credentials. The ZIP contains structures and a manifest with input IDs,
public source URLs, and SHA-256 hashes; a matching .manifest.json is written
beside it. The runtime verifies the structure hashes and preserves this
provenance in the per-input JSON. Use --help for all options and --force only
when intentionally replacing an existing prepared ZIP.

## Results and visualization

The **Summary** table has one row per input and five columns: structure, sites,
top probability, status, and structure file. Sites is the total qualifying count;
Truncated means the per-input site limit reduced the returned details.
The structure link opens the
first-model structure with predicted-site markers. The **Predicted sites** table
contains the site rank, probability, coordinating residues, template ID, and
template RMSD. Detailed coordinates and provenance remain available in the JSON
and complete results ZIP.

The Files panel has one canonical viewable structure and one PyMOL session per
input, with friendly names such as 01_1KLS_ZincSight.pdb and
01_1KLS_ZincSight.pse. Full mmCIF copies, per-input tables, and marker metadata
are collected inside the ZIP to keep the main result list compact. Tamarind uses
its own structure viewer; the session's PyMOL colors and labels are available
in the downloaded .pse.

| Output | Contents |
| --- | --- |
| 00_Summary.csv | Five-column summary with one row and viewable structure per input, including zero-site inputs |
| 01_Predicted_sites.csv | Six-column combined site table with probabilities, residues, and template matches |
| NN_label_ZincSight.pdb or .cif | One canonical structure per input: original first-model atoms and ions plus predicted-site markers |
| NN_label_ZincSight.pse | A ready-to-open PyMOL session per input |
| prediction.json | Full batch predictions, settings, provenance, resources, timings, and public/bundled artifact mappings |
| CITATION.txt | The ZincSight paper citation |
| ZincSight_results.zip | Complete detailed results, including full-precision mmCIF, PDB, PyMOL, per-input CSV/JSON, marker annotations, citation/BibTeX, and a public-to-bundle manifest |

Inside the ZIP, the detailed layout remains structures.csv, sites.csv,
sites/*.csv, predictions/*.json, and structures/*.{pdb,cif,pse} with matching
structures/*_sites.json annotations. The bundle also retains the complete
prediction.json and citation files. Public JSON distinguishes current
downloadable paths from bundle member paths; files described as bundle members
are opened after extracting the ZIP.

Coordinates remain in the input frame. The predicted Zn atoms are **location
markers**, not metal-identity assignments, and use a separate unused chain.
Original ions remain present. mmCIF and JSON retain full predicted coordinates;
standard PDB rounds coordinates to 0.001 Å. The public viewer uses PDB when it
can represent the input, and mmCIF otherwise. Long chain IDs or residue numbers
may need explicit mappings for PDB. Full mmCIF is always included in the bundle;
if the input exceeds PDB format limits, the result records the reason and uses
mmCIF/PyMOL.

PyMOL sessions show the original protein as a lavender cartoon, coordinating
residues as purple sticks with element-colored donor atoms, and predicted sites
as spheres colored from red to green by probability, labeled with site rank and
probability. Original atoms, experimental-metal selections, predicted markers,
and per-site coordinating selections remain distinguishable. Rotamer sampling
changes the prediction search; the session retains the original input
coordinates. A zero-site result still includes its structure and session, with
no invented markers and a header-only per-input sites CSV inside the bundle.

## Settings and compute

| Setting | Default and behavior |
| --- | --- |
| histidine_rotamers | true; samples alternative histidine orientations |
| minimum_probability | 0.0001; output filter matching the notebook, not a biological confidence threshold |
| max_sites | 100 per input; returns highest-probability sites, records the qualifying total and truncation |
| cpu_workers | 4; limited to the allocated CPUs and a 3 GiB budget per worker |

The deployment requests **4 vCPUs and 12 GiB memory**, with no GPU. The
whole batch shares one private embedded PostgreSQL (PGlite) instance and one
native ZincSight run. Available workers are capped by declared resources and
detected CPU/memory limits. This budgets 3 GiB per CPU worker; there is no GPU
implementation to accelerate this engine. Larger allocations can help native
parallel stages and batches, but speed depends on structure size and workload.
The platform job limit is 1,800 seconds; the scientific worker has a bounded
1,200-second timeout. Reduce the batch if it exceeds this limit.

## Validation and scientific status

The scientific engine and weights remain pinned to commit
7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03, with the source archive SHA-256
verified during the Docker build. The image applies a checked transport-only
patch so the job-local database uses a private Unix socket rather than loopback
TCP. Results record the engine revision, integration version 0.2.1, and
database_transport: unix_socket. To change the scientific revision, update the
Dockerfile pin/digest and ENGINE_REVISION together and rerun validation.

The **Tamarind validation** workflow validates the live Tamarind schema, runs
stdlib-only contract tests on the host, builds the dependency image, and runs
all Tamarind unit tests inside that image. Its integration check mounts the
entire integrations/tamarind directory and inputs over /app with networking
disabled. It retains exact agreement with the pinned wrapper for PDB/mmCIF,
rotamers on/off, truncation, and zero sites, then checks native file-list and
prepared-archive batches, duplicate basenames, the compact public file/table
contract, output ZIP integrity, citations, and reopening every bundled PyMOL
session. This still requires a separate live test
on Tamarind to verify the deployed platform's preview and download behavior.

~~~bash
# Host tests do not require BioPython or PyMOL.
for test_file in test_prepare_inputs.py test_tamarind.py test_input_bundle.py test_resources.py test_results_packager.py; do
  python3 -m unittest discover -s tests -p "$test_file" -v
done
docker build -t zincsight-tamarind:test .
python3 tests/check_tamarind_container.py \
  --image zincsight-tamarind:test --output /tmp/zincsight-tamarind-results
~~~

The inherited development preview contains reproducibility corrections that can
change historical predictions. Published accuracy benchmarks have not been
re-established for this revision. See [notebook regression](../mcp/NOTEBOOK_REGRESSION.md)
and [rotamer diagnosis](../mcp/ROTAMER_DIAGNOSIS.md). The adapter does not change
geometry/scoring formulas or trained weights. Predictions do not establish metal
identity, occupancy, or catalytic activity; absence of a prediction does not rule
out binding. The existing repository [license](../../LICENSE) applies.

