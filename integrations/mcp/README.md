# ZincSight for AI assistants (development preview)

**Historical-output equivalence is unresolved.** A rotamer-state fix now makes
all six notebook examples agree across native/embedded and single/batch runs.
Both the untouched original code and this revision still differ from the saved
October reference. See [the comparison and release blockers](NOTEBOOK_REGRESSION.md).

This optional MCP server exposes the existing ZincSight predictor to compatible
clients, including ChatGPT and Claude. It is **not yet a hosted service or a
published directory listing**. Geometry/scoring formulas and model weights are unchanged. The probability
step now returns cleanly for empty predictions and raises an explicit error for
scoring failures instead of continuing into invalid operations.
An optional embedded PostgreSQL backend supports task runtimes that cannot create
the operating-system user needed by pgserver. The standard backend remains the
default; two temporary-table helpers now explicitly drop their intermediate table.

**Reproducibility change:** motif IDs, coordinate/residue reads, template traversal
and greedy proximity grouping now use explicit SQL ordering. The previous queries
relied on unspecified database row order. A cross-backend test on 2CBA exposed
different candidates/probabilities, so this draft stabilizes that ordering in both
backends. It can change predictions compared with older versions; it needs broader
scientific regression review before merging. Backend comparisons below compare
the two modes of this revised pipeline, not all historical ZincSight outputs.
The histidine-rotamer helper also now keeps excluded indices local to each call;
previously one candidate could incorrectly remove conformers from a later one.

## Supported workflow

Ask: "Use ZincSight to predict possible zinc-binding sites in PDB 1A1F and report
the coordinating residues and site probabilities."

The `predict_zinc_binding_sites` tool accepts exactly one of:

- `pdb_id`: four-character accession; the server retrieves its PDB file from RCSB.
- `query_id`: a notebook-style PDB ID/assembly, UniProt accession, full AlphaFold
  model ID or MGYP identifier. PDB/AlphaFold inputs use mmCIF; ESM uses PDB.
  A UniProt accession resolves to F1 model v6, matching the notebook default.
- `structure_text`: actual PDB or mmCIF contents, with `file_format` set to `pdb`
  or `cif`. This can be an AlphaFold model. The assistant/client must have access
  to the actual contents; this preview does not implement a native upload widget.

Optional arguments: `histidine_rotamers` (default false, matching the CLI),
`minimum_probability` (default 0.0001, matching the public notebook and applied after prediction), `max_sites` (1–100,
default 50). The first structural model and all its chains are used. Limits:
5 MB input, 4000 protein residues, 180 seconds per prediction, one concurrent
prediction per server process. No sequence-to-structure prediction is performed.

Results include probabilities, zinc coordinates (angstrom), coordinating residue
identifiers, template IDs/RMSD, and a link to the software. Predictions do not
establish metal specificity, occupancy, catalytic function or experimental binding.
The default cutoff is the notebook's output filter, not a biological validation
threshold. Set `minimum_probability=0` (CLI: `--minimum-probability 0`) to inspect
all computed sites; the effective value is included in each result.

## Ubuntu setup in your own virtual environment

From the repository root, using a Python environment compatible with the existing
scientific requirements:

```bash
python3 -m venv .venv-mcp
source .venv-mcp/bin/activate
python -m pip install --upgrade pip
python -m pip install -r integrations/mcp/requirements.txt
python -m pip install pytest
python -m pytest tests/test_mcp.py -q
python -m integrations.mcp.server
```

The last command starts stdio MCP and waits for a client. Prediction dependencies
include PostgreSQL binaries supplied by pgserver and PyMOL imported by the existing
pipeline even when visualization is disabled. The legacy Biopython 1.78 pin may
need a C compiler on newer Python versions. Do not silently upgrade scientific
packages without regression checks against known ZincSight predictions.

A local MCP client can launch this command with its working directory set to the
repository root. For Claude Code, for example (replace both absolute paths):

```bash
claude mcp add zincsight -- /absolute/path/ZincSight/.venv-mcp/bin/python -m integrations.mcp.server
```

Run this from the repository root; clients started elsewhere should use a launcher
that changes into the repository first. No OpenAI or Anthropic API key is needed
by the predictor itself.

## Run inside a cloud task without permanent hosting

The experimental `pglite` backend uses PostgreSQL compiled to WebAssembly,
running inside the task via Node.js. It requires neither a Linux database user
nor a hosted database account. It does not bypass runtime permissions. The same
reference data, geometry formulas and probability model are used, with the explicit
ordering changes described above.
Install the Python requirements above, then use Node.js 22 or newer:

```bash
npm ci --prefix integrations/mcp/pglite --ignore-scripts --no-audit --no-fund
python -m integrations.mcp.predict --database pglite --pdb-id 1A1F
python -m integrations.mcp.predict --database pglite --query-id 2A0S-assembly1 --histidine-rotamers
python -m integrations.mcp.predict --database pglite --query-id P0A6G5 --histidine-rotamers
python -m integrations.mcp.predict --database pglite --structure protein.cif --format cif --histidine-rotamers
```

These commands print JSON and exit; no MCP client is needed. For an MCP client,
start the same server with this environment variable:

```bash
ZINCSIGHT_DATABASE_BACKEND=pglite python -m integrations.mcp.server
```

The CLI accepts `--minimum-probability` and returns up to 100 sites by default (`--max-sites` can change that). MCP
keeps its default of 50 and marks truncated output explicitly. mmCIF preparation
retains the original coordinate precision and selects the first model.

Each prediction creates a separate in-memory database and a loopback-only socket
on a dynamically allocated port. The worker, database and temporary files are
removed after the task, including the process group on timeout. Both Python and
Node must run in the same task environment.

PGlite has a single underlying PostgreSQL session. This adapter deliberately uses
autocommit for its disposable, one-structure job so nested connections cannot
hold each other waiting for a transaction. It explicitly releases the pipeline's
temporary motif table and loads the four bundled COPY blocks through PGlite's
blob API. It is not a general PostgreSQL server replacement or a shared batch
database. A failed prediction is discarded rather than returned as a valid result.

For one-skill installation, use the [prepared plugin package](DISTRIBUTION.md).
Its exported launcher has passed a fresh-install prediction test. The package
does not publish a directory listing or guarantee availability on every ChatGPT
plan/runtime; see the distribution guide for the remaining publication steps.

Embedded runtime tests (can run without a non-root OS account):

```bash
ZINCSIGHT_RUN_PGLITE=1 python -m pytest tests/test_mcp.py tests/test_mcp_integration.py tests/test_probability_empty_results.py -v
```

Enable `ZINCSIGHT_RUN_INTEGRATION=1` as well on a non-root Linux host to compare
native and embedded outputs on 1A1F, 2CBA and 1UBQ, with rotamers on and off.

PGlite documentation: https://pglite.dev/docs/about and
https://pglite.dev/docs/pglite-socket

## HTTP development and public launch

```bash
python -m integrations.mcp.server --transport streamable-http
```

The local endpoint is `http://127.0.0.1:8000/mcp`. For remote clients, deploy on a
persistent Linux host behind an HTTPS reverse proxy. Keep the application bound
to loopback. Configure the MCP SDK's allowed hosts/origins for your actual domain
rather than disabling DNS-rebinding protection. Apply proxy request-size limits,
rate limits and process/container CPU, memory and disk limits. Run as a non-root
service user. This preview has no built-in account authentication or durable job
queue; choose authentication and operating limits before public exposure.

Connect the HTTPS endpoint as a custom MCP server in ChatGPT, or as a custom
connector in Claude, subject to your account's available features. Test actual
client tool selection and returned output before requesting publication.

Each call receives its own temporary source copy, database, input and results.
The wrapper deletes these after the call and does not return internal logs.
RCSB sees accession download requests. Client chat history and infrastructure
logging have their own retention settings; review those before publishing a
privacy policy. Hard crashes may leave temporary files, so configure host cleanup.

### Public submission checklist

- Resolve historical scientific differences and extend validation beyond the
  six notebook examples; measure production memory use and request latency.
- For on-demand execution, package and fresh-install-test a skill and its isolated
  runtime. Deploy a persistent HTTPS endpoint only for the remote MCP route.
- Confirm code, model and template redistribution/service rights with the owners;
  no license is selected by this integration.
- Publish accurate support, privacy and terms pages with actual operator details.
- Prepare directory metadata and working example prompts; submit through the
  platform's current review process. Publication is subject to review, and does
  not guarantee recommendations to users who have not installed the plugin.

Suggested directory description:

> Predict potential zinc-binding sites from protein structures with ZincSight.
> Inspect coordinating residues, predicted zinc positions and probability scores.
> Accepts PDB accessions and PDB/mmCIF structure text, including AlphaFold models.

Do not advertise this preview as published or production-ready until these steps
are complete.

References:
- https://developers.openai.com/plugins/quickstart
- https://developers.openai.com/plugins/deploy/submission
- https://developers.openai.com/plugins/plugin-guidelines
- https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp

## Validation status (2026-10-30)

The [notebook-default comparison](NOTEBOOK_REGRESSION.md) now passes consistency
checks across backends and single/batch execution, but **not** equivalence to
the saved historical notebook outputs. Original-code reproductions also differ
from that reference. A configuration correction now applies the public
notebook's tracked 0.0001 probability cutoff. The earlier 137/138 totals were
unfiltered engine output; the corresponding public-input runs contain 73/70
sites after filtering. See the report for details.

**45 tests passed on Ubuntu 24.04 / Python 3.12 / Node.js 22 in GitHub Actions:**
https://github.com/MECHTI1/ZincSight/actions/runs/41036947949

- The standard entry point and both wrapper backends agree on PDB 1A1F, with
  histidine rotamer sampling both disabled and enabled.
- Native PostgreSQL and PGlite return identical site dictionaries for 1A1F,
  2CBA and 1UBQ in both rotamer modes, using this revision's explicit ordering.
- PDB and mmCIF representations return identical predictions for this control.
- An artificial negative control (all atom coordinates scaled tenfold) returns
  zero sites without crashing. This is a geometry control, not an experimentally
  established biological nonbinder.
- Unit tests cover empty-result handling and ensure a probability-calculation
  failure raises an error rather than reporting absence of binding.
- Wrapper checks cover input validation, parsing existing result CSVs, cleanup,
  timeouts, MCP metadata and an actual stdio MCP handshake/tool call.
- Embedded tests also cover simultaneous isolated positive/negative jobs and
  a real successful prediction through the MCP protocol.
- New regression tests cover rotamer exclusions leaking between calculations,
  notebook identifier resolution and preservation of mmCIF coordinate precision.
- An additional end-to-end workflow compares all six notebook examples using
  fresh public inputs. Native and embedded batches agree on all 137 returned
  sites, probabilities and positions. Calling the wrapper separately for each
  identifier, with rotamers enabled and no probability cutoff, gives the same
  results, including 8QEP and the UniProt-to-AlphaFold shorthand.

Embedded predictions also completed in the ChatGPT cloud development runtime,
where native pgserver could not create a Linux user. Python 3.12 and Node.js 24
were used there; both rotamer modes and PDB/mmCIF input were exercised.

These checks establish integration behavior on the tested cases, not new accuracy
benchmarks or equivalence to all older ZincSight outputs. The explicit-ordering
change above requires broader scientific regression review before merge. An
assistant skill and its exported plugin package now have runtime checks;
public-directory distribution, production memory limits and broader scientific
validation remain separate release work. Persistent
HTTPS hosting is necessary only if choosing the remote MCP deployment route.

Dependencies installed successfully with `CC=gcc` for legacy Biopython 1.78.
scikit-learn is pinned to 1.5.2 to match the bundled probability model. Run the
full suite on a non-root Linux host with RCSB access:

```bash
ZINCSIGHT_RUN_INTEGRATION=1 ZINCSIGHT_RUN_PGLITE=1 python -m pytest tests/test_mcp.py tests/test_mcp_integration.py tests/test_probability_empty_results.py tests/test_histidine_rotamer_state.py -v
```

Run the notebook consistency checks with fresh public structures:

```bash
python tests/notebook_regression.py --input-dir /tmp/zincsight-inputs --output-dir /tmp/zincsight-validation
```

The workflow retains unfiltered engine outputs for backend consistency checks,
but reference and default-wrapper comparisons apply the public notebook's
0.0001 cutoff. Use `--embedded-only` when native PostgreSQL is unavailable. An optional
`--expected /path/to/reference.csv` compares a local reference without uploading
it. CI publishes only newly computed results from public structure downloads.
