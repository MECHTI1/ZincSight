# ZincSight on Tamarind

**ZincSight predicts transition-metal binding sites in protein structures**, with
interpretable template matches, interactive structures and ready-to-open PyMOL
sessions. It was trained on zinc sites and benchmarked for zinc, iron, manganese,
cobalt, copper and nickel. The published study reports high precision and recall
for its transition-metal panel at probability ≥ 0.75; the separate zinc coverage
benchmark reached 97.7% sensitivity at probability ≥ 0.5.
See [the ZincSight paper](https://doi.org/10.1002/pro.70350) for the datasets and
benchmark conditions. Site markers locate binding sites without assigning the
bound metal identity.

## Run a batch

Upload many PDB/mmCIF files, one or more structure archives, or both. The tool
analyzes all protein chains in each input's first model and returns **every site
meeting the selected probability threshold**. There is no fixed structure-count,
per-structure byte, total-byte, protein-residue or returned-site cap.

For ID entry and file preparation, open
[Prepare Tamarind Batch in Colab](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/feature/tamarind-integration/integrations/tamarind/Prepare_Tamarind_Batch.ipynb).
Paste comma/space/newline-separated IDs, optionally upload ID lists or structures,
then prepare and download the input archives. Upload **all generated ZIP parts**
together under **Batch archives or prepared ID batches**. The notebook only
prepares inputs; predictions run on Tamarind. It embeds the preparer source and
needs no API key or GPU.

| Input | Accepted content |
| --- | --- |
| structures | One or many .pdb, .cif or .mmcif files |
| batch_archive | One or many ZIP, TAR, TAR.GZ or TGZ archives, including prepared ID ZIP parts; members may also be .ent or gzip-compressed structures |

Tamarind's upload transport accepts **200 MB per file**. The preparer splits
larger batches by compressed upload bytes, not by structure count, and includes
a hash-bearing manifest in every part. Each input is accounted for in the
aggregate manifest; nothing is silently omitted. A single structure whose
compressed archive cannot fit one upload produces a clear error. The actual
batch capacity is determined by the job's CPU, memory, disk and time allocation.

Supported IDs are PDB entries/assemblies, UniProt accessions, complete AlphaFold
model IDs, MGYP identifiers and complete TED domain IDs. UniProt resolves to
AlphaFold F1 model v6; explicit fragment/version IDs are preserved. Examples:
1KLS, 2A0S-assembly1, P12345, AF-P12345-F2-model_v4 and
AF-A0A287EUI1-F1-model_v4_TED01. ID availability depends on the public provider.
Tamarind jobs are offline, so download/ID resolution happens in the preparer.

The same standard-library helper works locally:

~~~bash
python3 integrations/tamarind/prepare_inputs.py \
  --ids "1KLS, 2A0S-assembly1, P12345" --output zincsight_inputs.zip

python3 integrations/tamarind/prepare_inputs.py \
  --ids-file ids.csv structures/ another_structure.pdb.gz \
  --output zincsight_inputs.zip
~~~

The CLI accepts .txt/.csv/.tsv ID lists, structures, folders and archives. It
spools structures to temporary disk, uses exact public-provider endpoint
allowlists, validates redirects and bounds download retries. Archive traversal,
links, special/encrypted members and path collisions are rejected. Optional
count, file-byte, total-byte, ID-list and archive/folder-entry guards are available
through --help and are **disabled by default**. --upload-max-bytes can choose a
smaller upload part size; --force replaces an existing prepared output family.
No account credentials or absolute local paths are included in the archives.

## Results

**Summary** has one row per input: structure, qualifying sites, top probability,
status and a structure link. **Predicted sites** lists all qualifying sites with
structure/source identity, site rank, probability, full X/Y/Z coordinates in Å,
coordinating residues, residue combination, template ID and template RMSD. The Files panel
keeps one friendly-named structure and PyMOL session per input, plus the complete
results ZIP and prediction JSON.

| Output | Contents |
| --- | --- |
| 00_Summary.csv | Five-column per-input summary, including zero-site structures |
| sites.csv | Complete site table: structure/source IDs, rank, probability, X/Y/Z coordinates, coordinating residues, residue combination and template matches |
| NN_label_ZincSight.pdb or .cif | One canonical structure per input, with original atoms/ions and predicted-site markers |
| NN_label_ZincSight.pse | A ready-to-open PyMOL session per input |
| prediction.json | Full predictions, settings, provenance, timings and public/bundle artifact mappings |
| CITATION.txt | ZincSight paper citation |
| ZincSight_results.zip | Full-precision mmCIF, PDB, PyMOL, per-input CSV/JSON, coordinates, marker annotations, mappings, citation/BibTeX and a public-to-bundle manifest |

The ZIP keeps the detailed structures.csv, sites.csv, sites/*.csv,
predictions/*.json and structures/*.{pdb,cif,pse} layout with matching marker
annotations. Public JSON distinguishes downloadable public paths from paths
inside the bundle. Coordinates stay in the input frame. PDB uses its standard
0.001 Å coordinate precision and may map long chain IDs/residue numbers;
full-precision mmCIF remains in the bundle. When PDB cannot represent the input,
the viewer uses mmCIF.

PyMOL shows a lavender protein cartoon, purple coordinating side chains with
element-colored donors, and predicted-site spheres colored by probability with
rank/probability labels. Original atoms, experimental metals and predictions
remain distinguishable. Rotamer sampling explores alternative histidine
orientations while visualization retains the original input coordinates.
Zero-site inputs still get a structure/session and no invented markers.

## Settings and resources

| Setting | Behavior |
| --- | --- |
| histidine_rotamers | Enabled by default |
| minimum_probability | Default 0.0001 matches BatchDrive; every qualifying site is returned |
| cpu_workers | Default 8; 1–8 workers within the available job allocation |

The tool requests **8 vCPUs, 24 GiB memory and 50 GiB working disk**, with a
**two-hour job budget**. The native CPU engine shares one private embedded
PostgreSQL/PGlite database across the batch. Workers are bounded by detected CPU
and memory availability, with 3 GiB budgeted per worker. Larger workloads may
need several jobs when they exceed the actual platform resources; failures are
reported rather than silently dropping structures or sites. The engine has no
GPU acceleration path.

## Deploy and validate

| Tamarind setting | Value |
| --- | --- |
| Repository | MECHTI1/ZincSight |
| Branch | feature/tamarind-integration |
| Source directory | Repository root |
| Existing/suggested tool ID | zincsight-preview |
| Display name | ZincSight |

Link the branch in [Tamarind Custom Tools](https://app.tamarind.bio/app/custom-tools),
build/deploy, test a batch with regular and zero-site inputs, then publish the
validated version to the organization. Verify the viewer, sessions and downloaded
bundle on Tamarind itself. Test/production jobs consume account compute.
See [Custom Tool conventions](https://docs.tamarind.bio/tamarind/custom-tools).

Integration revision 0.3.1 keeps the engine/weights pinned to
7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03 and its
source archive SHA-256 is checked during image build. Dependencies live under
/opt because Tamarind overlays the submitted source at /app. A checked
transport-only patch provides a job-private Unix database socket. Prediction and
packaging occur in private temporary storage; only final public files are copied
to /app/out, preserving Tamarind's logs and service-control files.

The validation workflow checks the live schema, notebook/helper equality,
standard-library contracts and image-dependent structure tests. Container checks
preserve exact agreement with the pinned wrapper, cover native and prepared
batches, and reopen every PyMOL session. A protected /app/out mount and open
platform-style log detect output-directory/service-file handling regressions.

~~~bash
for test_file in test_prepare_inputs.py test_tamarind.py test_input_bundle.py test_resources.py test_results_packager.py; do
  python3 -m unittest discover -s tests -p "$test_file" -v
done
docker build -t zincsight-tamarind:test .
python3 tests/check_tamarind_container.py \
  --image zincsight-tamarind:test --output /tmp/zincsight-tamarind-results
~~~

Published benchmark numbers describe the paper's evaluation. This integration's
pinned engineering revision includes reproducibility corrections; see
[notebook regression](../mcp/NOTEBOOK_REGRESSION.md) and
[rotamer diagnosis](../mcp/ROTAMER_DIAGNOSIS.md) for technical validation context.

Please cite Mechtinger G, Axel G, Kolodny R, Ben-Tal N. *Interpretable prediction
of zinc ion location in proteins with ZincSight*. Protein Science.
2025;34(11):e70350. [DOI](https://doi.org/10.1002/pro.70350).
The Paper link, CITATION.txt and bundled BibTeX reference this study; Tamarind
controls its native Citations tab separately. The repository
[license](../../LICENSE) applies.
