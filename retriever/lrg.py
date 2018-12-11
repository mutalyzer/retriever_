import os
from urllib.error import URLError
from urllib.request import urlopen

# Prefix URL from where LRG files are fetched.
LRG_PREFIX_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/lrgex/'

# Maximum size for uploaded and downloaded files (in bytes).
MAX_FILE_SIZE = 10 * 1048576  # 10 MB


def fetch_lrg(reference_id, size_on=True):
    """
    Fetch the LRG file content.

    :param size_on: flag for the maximum sequence length
    :param reference_id: the name of the LRG file to fetch
    :returns: the file content or None in case of an error
    """
    url = '{}/{}.xml'.format(LRG_PREFIX_URL, reference_id)
    filename = None

    try:
        handle = urlopen(url)
    except URLError:
        return None

    info = handle.info()

    if (info['Content-Type'] == 'application/xml' and
            'Content-length' in info):
        if size_on:
            length = int(info['Content-Length'])
            if 512 > length  or length > MAX_FILE_SIZE:
                print('4, EFILESIZE, Filesize {} is not within the allowed '
                      'boundaries (512 < filesize < {} ) for {}.'
                      .format(length, MAX_FILE_SIZE // 1048576, reference_id))
                handle.close()
                return None
        raw_data = handle.read()
        handle.close()
        return raw_data
    else:
        print('4, ERECPARSE, This is not an LRG record.')
        handle.close()
        return None
