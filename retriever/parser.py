from .parsers import genbank, lrg


def get_reference_type(content):
    reference_type = None

    return reference_type


def parse(content, reference_type=None):

    print(content[0:100], reference_type)

    if reference_type is None:
        reference_type = get_reference_type(content)

    if reference_type == 'lrg':
        reference = lrg.parse(content)
    elif reference_type == 'genbank_ncbi':
        reference = genbank.parse(content)
    else:
        return None

    return reference.loci_to_json_model()
