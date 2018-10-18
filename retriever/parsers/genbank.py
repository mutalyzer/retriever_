import hashlib
import io
import os

from Bio import SeqIO, SeqRecord, SeqFeature
from ..reference import Position, Locus, Reference

FEATURES = {
    'source': {
        'qualifiers': [
            'organism',
            'mol_type',
            'chromosome',
            'organelle'
        ]
    },
    'gene': {
        'qualifiers': [
            'gene',
            'gene_synonym',
            'db_xref'
        ]
    },
    'mRNA': {
        'qualifiers': [
            'transcript_id',
            'locus_tag',
            'product',
            'gene',
            'gene_synonym'
        ]
    },
    **(dict.fromkeys(
        ['mRNA', 'misc_RNA', 'ncRNA', 'tRNA', 'rRNA', 'precursor_RNA'], {
            'qualifiers': [
                'transcript_id',
                'locus_tag',
                'product',
                'gene',
                'gene_synonym'
            ]
        }
    )),
    'CDS': {
        'qualifiers': [
            'gene',
            'protein_id',
            'codon_start',
            'transl_table',
            'product',
            'db_xref'
        ]
    }
}


def parse(content):

    record = SeqIO.read(io.StringIO(content), 'genbank')

    # The provided record must be an instance of BioPython SeqRecord.
    if not isinstance(record, SeqRecord.SeqRecord):
        print('Record_input not of BioPython SeqRecord type.')
        return None

    # If there are no features in the record we can stop.
    if not record.features:
        return None

    reference = Reference()

    extract_features(reference, record)

    print(reference)

    return reference


def extract_features(reference, record):
    """
    Loops over the record features and extracts them into Locus instances.

    :arg Reference reference: The actual reference.
    :arg SeqRecord record: A biopython record.
    """
    for feature in record.features:
        if feature.type in FEATURES:
            configuration = FEATURES[feature.type]
            locus = convert_biopython_feature(feature, configuration)
            reference.loci.append(locus)
    print(reference)


def convert_biopython_feature(biopython_feature, config=None):
    """
    Gathers all the relevant information from a BioPython feature instance
    and populates the Locus attributes based on the provided
    configuration.

    It checks also if the feature is composed of multiple parts and if so
    adds them to the sequence 'parts' attribute list.

    :param biopython_feature: The BioPython feature from where the Locus
                              information is to be extracted.
    """
    # TODO: Could be converted into a static method or a separate package.
    if not isinstance(biopython_feature, SeqFeature.SeqFeature):
        return
    if config is None:
        print('No config.')
        return
    if not isinstance(config, dict):
        print('Config not a dictionary.')
        return
    if 'qualifiers' not in config:
        print('No qualifiers.')
        return

    start = Position(str(biopython_feature.location.start))
    start.add_int(1)
    end = Position(str(biopython_feature.location.end))
    locus_type = biopython_feature.type
    orientation = biopython_feature.strand

    qualifiers = {}
    # Get all the qualifiers.
    for qualifier in config['qualifiers']:
        qualifiers[qualifier] = extract_qualifier(qualifier, biopython_feature)

    parts = None
    # Check if it is a compound sequence.
    if len(biopython_feature.location.parts) > 1:
        parts = []
        # Gather all the parts. The part type cannot be gathered now as there
        # is no information on that. Should be determined later.
        for part in biopython_feature.location.parts:
            part_start = Position(str(part.start))
            part_start.add_int(1)
            part_sequence = Locus(start=part_start,
                                  end=str(Position(part.end)),
                                  orientation=orientation)
            parts.append(part_sequence)

    locus = Locus(start=start, end=end, orientation=orientation,
                  locus_type=locus_type, parts=parts)
    return locus


def extract_qualifier(qualifier, biopython_feature):
    """
    Extracts a specific qualifier from a biopython feature and adds it to
    the sequence qualifiers.
    If there is no such qualifier none is returned.

    :param qualifier: the qualifier to be extracted.
    :param biopython_feature: the feature to which the qualifier belongs to.
    :return: True if qualifier present in the feature and False otherwise.
    """
    if qualifier in biopython_feature.qualifiers:
        return biopython_feature.qualifiers[qualifier][0]


def extract_dbxref_qualifiers(qualifier, parts, biopython_feature):
    """
    Extracts qualifiers which can be part of the same biopython feature, as
    for example the "dbxref" which can be as follows:
    - /db_xref="HGNC:HGNC:26049"
    - /db_xref="CCDS:CCDS53487.1"
    - /db_xref="GeneID:55657"

    :param qualifier: The qualifier for which the extraction is performed.
    :param parts: What parts to extract.
    :param biopython_feature: The biopython feature from which the
                              qualifier is extracted.
    """
    qualifiers = {}
    for part in parts:
        qualifiers[qualifier + '_' + part] = None
    if qualifier in biopython_feature.qualifiers:
        for part in parts:
            for ref in biopython_feature.qualifiers[qualifier]:
                if part == ref.split(':')[0]:
                    qualifiers[qualifier + '_' + part] = ':'.\
                        join(ref.split(':')[1:])
