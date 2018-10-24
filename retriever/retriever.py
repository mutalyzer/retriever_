from .lrg import fetch_lrg
from .ncbi import fetch_ncbi
from .parsers import genbank

from pathlib import Path

CACHE_PATH = '/tmp'
CACHE = True


def retrieve(reference_id, size_on=True, parse=False):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods (lrg or ncbi genbank so far).

    :arg bool parse: Flag for parsing or not the reference.
    :arg str reference_id: The id of the reference.
    :arg bool size_on: Flag for the maximum sequence length.
    :return: The reference raw content or its equivalent parse tree serialized.
    :rtype: str
    """
    content = None
    if CACHE:
        path = Path(CACHE_PATH) / reference_id
        if path.is_file():
            with path.open() as f:
                content = f.read()
    if content is None:
        if 'LRG' in reference_id:
            content = fetch_lrg(reference_id, size_on)
        else:
            content = fetch_ncbi(reference_id, size_on)

    if parse:
        reference = genbank.parse(content)
        if reference and CACHE:
            path = Path(CACHE_PATH) / reference_id
            with path.open('w') as f: f.write(content)
        return reference

