from .lrg import fetch_lrg
from .ncbi import fetch_ncbi
from .parsers import genbank


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

    if 'LRG' in reference_id:
        content = fetch_lrg(reference_id, size_on)
    else:
        content = fetch_ncbi(reference_id, size_on)

    if parse:
        return genbank.parse(content)
