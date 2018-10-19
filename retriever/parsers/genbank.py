import hashlib
import io
import os

from Bio import SeqIO, SeqRecord, SeqFeature
from ..reference import Position, Locus, Reference

KEY_FEATURES = {
    'gene': {
        'qualifiers': {
            'gene',
            'gene_synonym',
            'db_xref'
        },
        'key': 'gene'
    },
    'mRNA': {
        'qualifiers': {
            'transcript_id',
            'locus_tag',
            'product',
            'gene',
            'gene_synonym'
        },
        'key': 'transcript_id'
    },
    **(dict.fromkeys(
        ['mRNA', 'misc_RNA', 'ncRNA', 'tRNA', 'rRNA', 'precursor_RNA'], {
            'qualifiers': {
                'transcript_id',
                'locus_tag',
                'product',
                'gene',
                'gene_synonym'
            },
        'key': 'transcript_id'
        }
    )),
    'CDS': {
        'qualifiers': {
            'gene',
            'protein_id',
            'codon_start',
            'transl_table',
            'product',
            'db_xref'
        },
        'key': 'protein_id'
    }
}

SOURCE_FEATURES = {
    'qualifiers': {
        'organism',
        'mol_type',
        'chromosome',
        'organelle'
    }
}


def parse(content):
    """
    Parses a genbank file and returns the defined model.

    :arg str content: Genbank file content.
    :return: The corresponding reference instance.
    :rtype: Reference
    """
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

    return reference


def extract_features(reference, record):
    """
    Loops over the record features and extracts them into Locus instances.

    :arg Reference reference: The actual reference.
    :arg SeqRecord record: A biopython record.
    """
    for feature in record.features:
        if feature.type == 'source':
            reference.source = convert_biopython_feature(feature,
                                                         SOURCE_FEATURES)
        elif feature.type in KEY_FEATURES:
            config = KEY_FEATURES[feature.type]
            locus = convert_biopython_feature(feature, config)
            key = config.get('key')
            reference.loci[locus.qualifiers[key]] = locus


def convert_biopython_feature(biopython_feature, config=None):
    """
    Gathers all the relevant information from a BioPython feature instance
    and populates the Locus attributes based on the provided
    configuration.

    It checks also if the feature is composed of multiple parts and if so
    adds them to the sequence 'parts' attribute list.

    :arg SeqFeature biopython_feature: The BioPython feature from where the
    Locus information is to be extracted.
    :arg dict config: Configuration dictionary for the extraction process.
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
        print('No qualifiers in configuration.')
        return

    start = Position(str(biopython_feature.location.start))
    start.add_int(1)
    end = Position(str(biopython_feature.location.end))
    locus_type = biopython_feature.type
    orientation = biopython_feature.strand

    qualifiers = extract_qualifiers(biopython_feature.qualifiers, config)

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
                  locus_type=locus_type, parts=parts, qualifiers=qualifiers)
    return locus


def extract_qualifiers(biopython_qualifiers, config):
    """
    Extracts the biopython qualifiers mentioned in the configuration.

    :param biopython_qualifiers: The qualifiers to be extracted.
    :param config: The corresponding configuration based on which the
    extraction is performed.
    :return: Extracted qualifiers.
    :rtype: dict
    """
    qualifiers = {}
    for k in biopython_qualifiers.keys():
        if k in config['qualifiers']:
            if k == 'db_xref':
                qualifiers.update(
                    extract_dbxref_qualifiers(biopython_qualifiers[k]))
            else:
                qualifiers[k] = biopython_qualifiers[k][0]
    return qualifiers


def extract_dbxref_qualifiers(biopython_qualifier):
    """
    Extracts qualifiers which can be part of the same biopython feature, as
    for example the "dbxref" which can be as follows:
    - /db_xref="HGNC:HGNC:26049"
    - /db_xref="CCDS:CCDS53487.1"
    - /db_xref="GeneID:55657"

    :arg str biopython_qualifier: The qualifier for which the extraction is
    performed.
    :return: Dictionary containing the multiple db_xref fields.
    :rtype: dict
    """
    qualifiers = {}
    for sub_qualifier in biopython_qualifier:
        if 'HGNC' in sub_qualifier:
            qualifiers['HGNC'] = sub_qualifier.rsplit(':', 1)[1]
        elif 'CCDS' in sub_qualifier:
            qualifiers['CCDS'] = sub_qualifier.rsplit(':', 1)[1]
        elif 'GeneID' in sub_qualifier:
            qualifiers['GeneID'] = sub_qualifier.rsplit(':', 1)[1]
    return qualifiers
