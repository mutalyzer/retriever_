from .lrg import fetch_lrg
from .ncbi import fetch_ncbi, get_reference_summary
from .parsers import genbank, lrg

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
    checks = {
        'ncbi': False,
        'lrg': False
    }

    content = None
    if CACHE:
        path = Path(CACHE_PATH) / reference_id
        if path.is_file():
            with path.open() as f:
                content = f.read()

    if 'LRG' in reference_id:
        content = fetch_lrg(reference_id, size_on)
        checks['lrg'] = True
        reference_type = 'lrg' if content else None
    else:
        content, check = fetch_ncbi(reference_id, size_on)
        if check:
            checks['ncbi'] = True
        reference_type = 'genbank_ncbi' if content else None

    if (content is None) and (not checks['lrg']):
        content = fetch_lrg(reference_id, size_on)
        checks['lrg'] = True
    if (content is None) and (not checks['ncbi']):
        content, check = fetch_ncbi(reference_id, size_on)

    print(content)

    if parse:
        if reference_type == 'lrg':
            reference = lrg.parse(content)
        elif reference_type == 'genbank':
            reference = genbank.parse(content)
            if reference and CACHE:
                path = Path(CACHE_PATH) / reference_id
                with path.open('w') as f:
                    f.write(content)

        reference.loci_to_json_model()

        return reference
