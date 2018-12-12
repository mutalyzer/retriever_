from .lrg import fetch_lrg
from .ncbi import fetch_ncbi, NcbiConnectionError

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
    :return: The reference raw content and its type.
    :rtype: tuple
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
            # Todo: Retrieve the reference type also
            return content, None

    if 'LRG' in reference_id:
        content = fetch_lrg(reference_id, size_on)
        checks['lrg'] = True
        reference_type = 'lrg' if content else None
    else:
        try:
            content = fetch_ncbi(reference_id, size_on)
        except NcbiConnectionError:
            print('NCBI connection error.')
        else:
            checks['ncbi'] = True
            reference_type = 'genbank_ncbi' if content else None

    if (content is None) and (not checks['lrg']):
        content = fetch_lrg(reference_id, size_on)
        reference_type = 'lrg'
    if (content is None) and (not checks['ncbi']):
        try:
            content = fetch_ncbi(reference_id, size_on)
        except NcbiConnectionError:
            print('NCBI connection error.')
            return content, None
        else:
            checks['ncbi'] = True
            reference_type = 'genbank_ncbi' if content else None

    return content, reference_type
