from http.client import HTTPException
from urllib.error import HTTPError

from Bio import Entrez

# Maximum size for uploaded and downloaded files (in bytes).
MAX_FILE_SIZE = 10 * 1048576  # 10 MB

# This email address is used in contact information on the website and sent
# with NCBI Entrez calls.
EMAIL = 'mutalyzer@lumc.nl'

Entrez.email = EMAIL


def length(reference_id):
    """
    Retrieves a the sequence length.

    :param size_on: flag for the maximum sequence length
    :param reference_id: the id of the reference
    :returns: sequence length
    """
    try:
        net_handle = Entrez.esummary(db='nuccore', id=reference_id)
        record = Entrez.read(net_handle)
        net_handle.close()
    except (IOError, HTTPError, HTTPException) as e:
        print('-1, INFO, Error connecting to Entrez nuccore database: {}.'
              .format(e))
        print('4, ERETR, Could not retrieve record length for {}.'
              .format(reference_id))
        return None

    return int(record[0]['Length'])


def fetch_ncbi(reference_id, size_on=True):
    """
    Retrieves a genbank reference from the NCBI.

    :param reference_id: the id of the reference
    :param size_on: consider or not the maximum sequence length
    :returns: reference content
    """
    if size_on:
        if length(reference_id) is None:
            print('Retrieval stopped since there is no record length.')
            return
        if length(reference_id) > MAX_FILE_SIZE:
            print('4, ERETR, Could not retrieve {} (exceeds maximum file '
                  'size of {} megabytes).'
                  .format(reference_id, MAX_FILE_SIZE // 1048576))
            return None
    try:
        net_handle = Entrez.efetch(
            db='nuccore', id=reference_id, rettype='gbwithparts',
            retmode='text')
        raw_data = net_handle.read()
        net_handle.close()
    except (IOError, HTTPError, HTTPException) as e:
        print('-1, INFO, Error connecting to Entrez nuccore database: {}.'
              .format(e))
        print('4, ERETR, Could not retrieve record {}.'
              .format(reference_id))
        return None

    return raw_data
