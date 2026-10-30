# ZincSight for AI assistants (development preview)

This optional MCP server exposes the existing ZincSight predictor to compatible
clients, including ChatGPT and Claude. It is **not yet a hosted service or a
published directory listing**. The scientific predictor is unchanged.

## Supported workflow

Ask: "Use ZincSight to predict possible zinc-binding sites in PDB 1A1F and report
the coordinating residues and site probabilities."

The `predict_zinc_binding_sites` tool accepts exactly one of:

- `pdb_id`: four-character accession; the server retrieves its PDB file from RCSB.
- `structure_text`: actual PDB or mmCIF contents, with `file_format` set to `pdb`
  or `cif`. This can be an AlphaFold model. The assistant/client must have access
  to the actual contents; this preview does not implement a native upload widget.

Optional arguments: `histidine_rotamers` (default false, matching the CLI),
`minimum_probability` (default 0, applied after prediction), `max_sites` (1–100,
default 50). The first structural model and all its chains are used. Limits:
2 MB input, 2000 protein residues, 180 seconds per prediction, one concurrent
prediction per server process. No sequence-to-structure prediction is performed.

Results include probabilities, zinc coordinates (angstrom), coordinating residue
identifiers, template IDs/RMSD, and a link to the software. Predictions do not
establish metal specificity, occupancy, catalytic function or experimental binding.
There is no universally selected probability cutoff in this wrapper.

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

- Validate real positive, negative, rotamer and mmCIF cases against the original
  predictor; confirm output equivalence and measure cold/warm request latency.
- Deploy and test the persistent HTTPS endpoint, concurrency and failure behavior.
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

The wrapper tests cover validation, parsing the repository's existing prediction
CSV, failure-versus-no-sites handling, temporary-directory cleanup, timeout cleanup,
tool metadata, and an actual stdio MCP client handshake/tool call. These are not
scientific accuracy or end-to-end prediction tests.

Installation succeeded with Python 3.12 using `CC=gcc` to build the legacy
Biopython dependency. The integration pins scikit-learn 1.5.2 to match the version
recorded in the bundled probability model. Importing the prediction pipeline
succeeded. A real 1A1F prediction could not complete in the development environment:
pgserver attempted `useradd ... pgserver` under root and the environment rejected
that operation. Validate on a normal non-root Linux host before launch.
