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

Under **Protein structures or batch archives**, upload many PDB/mmCIF files,
one or more structure archives, or both together in the same picker. The tool
analyzes all protein chains in each input's first model and returns **every site
meeting the selected probability threshold**. There is no fixed structure-count,
per-structure byte, total-byte, protein-residue or returned-site cap.

### Fetch IDs inside Tamarind — no Colab

Tamarind's native [Database](https://app.tamarind.bio/database) can fetch
coordinate-bearing structures before ZincSight's offline container starts:

1. Open **Database → +New → Fetch by identifier**.
2. Choose **RCSB** for experimental PDB IDs such as `1KLS`, or **AlphaFold**
   for UniProt accessions such as `P69905` and model entries such as
   `AF-A0A3W0S055-F1`. Paste one or many IDs separated by spaces, commas or
   newlines. Choose the source explicitly; add each source's IDs separately.
3. Review the fetched structure previews and import them into a named group.
4. In that group, use **Export Structures** to download a ZIP, then unzip it.
5. Open ZincSight and multi-select the exported **PDB files** under
   **Protein structures or batch archives**.

This handoff was verified with all three example structures: the PDB uploads
completed, Submit became enabled, and the native 3D preview loaded `1KLS`.
No new prediction was submitted. This is a native Database workflow followed by
file upload; it does not add inline ID/CSV/Excel tabs to ZincSight.

For this native-export route, unzip before choosing the files. Uploading its ZIP
directly in the tested chooser opened **Import PDB Files from ZIP**, which asks
to create another molecule group; that dialog was canceled. Existing ordinary
structure/archive inputs remain supported by the adapter.

The exported files preserved every first-model atom and its coordinates in the
three checked examples, including the original `1KLS` Zn ion, with no residue
cropping. Export uses PDB serialization and model 1; `1KLS` models 2–30 are not
retained. ZincSight also analyzes model 1.

Native AlphaFold entry IDs differ from filename-style IDs such as
`AF-A0A3W0S055-F1-model_v6`. The native exporter omits model-version headers;
its checked coordinates match the v6 references, but that is not independent
proof of the fetched release. For an explicitly required fragment/version, use
the version-preserving reader below or upload the original provider coordinate
file. See [Database documentation](https://docs.tamarind.bio/tamarind/database).

### Prepare ID-list files locally — optional

**TXT, CSV and TSV structure-ID lists are supported by the input preparer.**
It downloads the requested structures before the offline Tamarind prediction job.
Run the standard-library Python helper below on a computer with internet access;
Colab is optional. Upload **all generated ZIP parts** together under **Protein
structures or batch archives**. An ID-only text file does not contain the
coordinates required by the Tamarind runtime and is not a direct prediction input.

| Input | Accepted content |
| --- | --- |
| structures | One or many .pdb, .cif, .mmcif or .ent structures, ZIP/TAR/TAR.GZ/TGZ archives, gzip-compressed structures, or any mixture; includes prepared ID ZIP parts |

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
For this optional reader, download/ID resolution happens locally before the offline job; the native Database route above is a separate option.

Run an ID list locally, without Colab or an API key:

~~~bash
python3 integrations/tamarind/prepare_inputs.py \
  --ids-file my_structures.tsv --output zincsight_inputs.zip
~~~

Plain TXT lists may contain commas, spaces, tabs or newlines, including quoted
IDs. For example:

~~~text
1KLS, P69905, "AF-A0A3W0S055-F1-model_v6.cif"
2A0S-assembly1
~~~

CSV/TSV tables may use an ID column such as `structure_id`, `model_id`,
`alphafold_model_id`, `uniprot_accession` or UniProt's `Entry`; descriptive columns
are ignored. The first recognized ID column is used. Without a recognized header,
every nonempty cell is treated as ID text. Separators are detected from the content,
so a comma-separated list saved as `.tsv` also works. `#` introduces a comment.
See the preserved [example lists](https://github.com/MECHTI1/ZincSight/tree/c0644e3798465d99ee92b5357537645e2f075e58/integrations/tamarind/examples) for ready-to-edit TXT, CSV and TSV files.

Known structure filename suffixes (`.pdb`, `.cif`, `.mmcif`, `.ent`, optionally
followed by `.gz`) are accepted as ID spellings. For example,
`AF-P12345-F2-model_v4.cif.gz` still requests fragment **F2**, version **v4**;
the reader does not replace an explicitly requested model with F1/v6. URLs and
filesystem paths are not structure IDs; pass local files separately.

The helper also accepts pasted IDs and mixtures of ID lists and local structures:

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

For an optional browser notebook, use
[Prepare Tamarind Batch in Colab](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/c0644e3798465d99ee92b5357537645e2f075e58/integrations/tamarind/Prepare_Tamarind_Batch.ipynb).
It embeds the same helper and accepts the same ID lists, structures and archives.
It only prepares inputs; predictions run on Tamarind. No API key or GPU is needed.

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
| NN_label_ZincSight_structure.pdb or .cif | One canonical structure per input, with original atoms/ions and predicted-site markers |
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

The single declared structure output matches `*_ZincSight_structure.*`, covering
both PDB and mmCIF while keeping input copies and separately named PyMOL sessions
outside the structure output contract.

PyMOL shows a lavender protein cartoon, purple coordinating side chains with
element-colored donors, and predicted-site spheres colored by probability with
rank/probability labels. The default **Overview** scene shows the highest-scoring
site's label while retaining every predicted sphere. **All_site_labels** shows
all stored site labels; switching back to **Overview** restores the compact view.
Original atoms, experimental metals and predictions
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
| Branch | master |
| Source directory | Repository root |
| Existing tool ID | transition-metal-binding-site-predictor |
| Display name | ZincSight |

Use the existing university tool in [Tamarind Custom Tools](https://app.tamarind.bio/app/custom-tools),
with its Git source tracking `master`; do not recreate it. Build/deploy,
test a batch with regular and zero-site inputs, then publish the
validated version to the organization. Verify the viewer, sessions and downloaded
bundle on Tamarind itself. Test/production jobs consume account compute.
See [Custom Tool conventions](https://docs.tamarind.bio/tamarind/custom-tools).

Integration revision 0.3.4 keeps the engine/weights pinned to
7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03 and its
source archive SHA-256 is checked during image build. Dependencies live under
/opt because Tamarind overlays the submitted source at /app. A checked
transport-only patch provides a job-private Unix database socket. Prediction and
packaging occur in private temporary storage; only final public files are copied
to /app/out, preserving Tamarind's logs and service-control files.

The tested engine and integration revisions remain ancestors of `master`.
The image downloads the immutable engine commit, including its database launcher
and dependency lock; those preview files do not need to be copied into the
current `master` tree. Keep the engine revision and archive digest unchanged
unless a separately validated engine update is intended. The optional preparation
notebook and example files above are also retained at immutable commit links.

The validation workflow runs on `master` and checks the live schema, adapter syntax,
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
[notebook regression](https://github.com/MECHTI1/ZincSight/blob/7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03/integrations/mcp/NOTEBOOK_REGRESSION.md) and
[rotamer diagnosis](https://github.com/MECHTI1/ZincSight/blob/7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03/integrations/mcp/ROTAMER_DIAGNOSIS.md) for technical validation context.

Please cite Mechtinger G, Axel G, Kolodny R, Ben-Tal N. *Interpretable prediction
of zinc ion location in proteins with ZincSight*. Protein Science.
2025;34(11):e70350. [DOI](https://doi.org/10.1002/pro.70350).
The Paper link, CITATION.txt and bundled BibTeX reference this study; Tamarind
controls its native Citations tab separately. The repository
[license](../../LICENSE) applies.
