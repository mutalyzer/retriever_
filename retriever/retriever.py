from .lrg import fetch_lrg
from .ncbi import fetch_ncbi


def retrieve(reference_id, size_on=True):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods (lrg or ncbi genbank so far).

    :arg str reference_id: The id of the reference.
    :arg bool size_on: Flag for the maximum sequence length.
    :return: The actual reference content in the corresponding format.
    :rtype: str
    """
    if 'LRG' in reference_id:
        return fetch_lrg(reference_id, size_on)
    else:
        return fetch_ncbi(reference_id, size_on)