from .parsers import genbank, lrg


def get_reference_type(content):
    if content.startswith('<?xml version='):
        return 'lrg'
    elif content.startswith('LOCUS'):
        return 'genbank_ncbi'


def parse(content, reference_type=None):

    if reference_type is None:
        reference_type = get_reference_type(content)

    if reference_type == 'lrg':
        reference = lrg.parse(content)
    elif reference_type == 'genbank_ncbi':
        reference = genbank.parse(content)
    else:
        return None

    return reference
