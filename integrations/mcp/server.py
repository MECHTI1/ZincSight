"""Run with python -m integrations.mcp.server [--transport streamable-http]."""
import argparse
import threading
from typing import Literal

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from integrations.mcp.runner import run_prediction

mcp = FastMCP("ZincSight", host="127.0.0.1", port=8000, stateless_http=True, json_response=True)
_capacity = threading.BoundedSemaphore(1)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                                     idempotentHint=True, openWorldHint=True))
def predict_zinc_binding_sites(
    pdb_id: str | None = None, structure_text: str | None = None,
    file_format: Literal["pdb", "cif"] = "pdb", histidine_rotamers: bool = False,
    minimum_probability: float = 0.0, max_sites: int = 50,
) -> dict:
    """Predict possible zinc-binding sites from an existing protein structure using ZincSight.

    Supply exactly one: a four-character PDB ID (retrieved from RCSB), or the actual
    PDB/mmCIF text (including AlphaFold models). Never invent structure coordinates.
    File paths, upload IDs, sequences and arbitrary URLs are not accepted.
    Maximum 2 MB / 2000 protein residues. Uses the first model, all chains.
    Optional histidine rotamer sampling; default off matches the existing CLI.
    Returns probabilities, coordinating residues, coordinates in angstrom, template
    references and limitations. minimum_probability filters returned sites only.
    This does not establish binding experimentally, metal specificity or catalysis.
    Inputs are processed in temporary server files, deleted after each call.
    """
    if not _capacity.acquire(blocking=False):
        raise ValueError("ZincSight is processing another request. Please retry shortly.")
    try:
        return run_prediction(pdb_id=pdb_id, structure_text=structure_text,
                              file_format=file_format, histidine_rotamers=histidine_rotamers,
                              minimum_probability=minimum_probability, max_sites=max_sites)
    finally:
        _capacity.release()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport", choices=["stdio", "streamable-http"], default="stdio")
    args = parser.parse_args()
    mcp.run(transport=args.transport)
