"""
Genbank file parser.

Creates a reference model based on genbank file.

Additionaly it should also provide some checks:
- what features do not contain the the configuration qualifiers.
- what features were not linked.
- what features do not have a key.
"""
import hashlib
import io
import os
import json

from Bio import SeqIO, SeqRecord, SeqFeature
from ..reference import Position, Locus, Reference
from ..ncbi import link_reference, compose_reference, decompose_reference

FEATURES = {
    'gene': {
        'qualifiers': {
            'gene',
            'gene_synonym',
            'db_xref'
        },
        'key': 'gene'
    },
    **(dict.fromkeys(
        ['mRNA', 'misc_RNA', 'ncRNA', 'tRNA', 'rRNA', 'precursor_RNA'], {
            'qualifiers': {
                'transcript_id',
                'locus_tag',
                'product',
                'gene',
                'gene_synonym',
                'db_xref'
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
            'db_xref',
            'coded_by'
        },
        'key': 'protein_id'
    },
    'Protein': {
        'qualifiers': {
            'product'
        }
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

ANNOTATIONS = {
    'accessions',
    'sequence_version',
    'molecule_type',
    'date'
}


def parse(content):
    """
    Parses a genbank string and returns the defined reference model.

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

    loci = _extract_features(reference, record)

    info = _build_info(record, reference, loci)

    if info.get('mol_type') == 'p':
        _add_extra_non_feature_loci(loci, info.get('id'))

    _construct_dependencies(loci)

    reference.loci = loci

    reference.info = info

    reference.type = 'genbank'

    return reference


def _extract_features(reference, record):
    """
    Loops over the record features and extracts them into Locus instances.

    :arg Reference reference: The actual reference.
    :arg SeqRecord record: A biopython record.
    """

    loci = {}

    for feature in record.features:
        if feature.type == 'source':
            reference.source = _convert_biopython_feature(feature,
                                                          SOURCE_FEATURES)
        elif feature.type in FEATURES:
            config = FEATURES[feature.type]
            locus = _convert_biopython_feature(feature, config)
            key = config.get('key')
            if key and locus.qualifiers.get(key):
                if feature.type not in loci:
                    loci[feature.type] = {}
                loci[feature.type][locus.qualifiers[key]] = locus
                reference.loci[locus.qualifiers[key]] = locus
            else:
                if loci.get('keyless') is None:
                    loci['keyless'] = {}
                if feature.type not in loci['keyless']:
                    loci['keyless'][feature.type] = []
                loci['keyless'][feature.type].append(locus)
    return loci


def _convert_biopython_feature(biopython_feature, config=None):
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

    qualifiers = _extract_qualifiers(biopython_feature.qualifiers, config)

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


def _extract_qualifiers(biopython_qualifiers, config):
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
                    _extract_dbxref_qualifiers(biopython_qualifiers[k]))
            else:
                qualifiers[k] = biopython_qualifiers[k][0]
    return qualifiers


def _extract_dbxref_qualifiers(biopython_qualifier):
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


def _extract_mol_type(reference, loci):
    """
    Extract molecule type for the reference.
    """
    mol_type = None
    if reference.source.qualifiers.get('mol_type'):
        if reference.source.qualifiers['mol_type'][0] in ['mRNA', 'transcribed RNA']:
            mol_type = 'n'
        elif reference.source.qualifiers['mol_type'][0] in ['genomic DNA']:
            mol_type = 'g'
    if reference.source.qualifiers.get('organelle'):
        if reference.source.qualifiers['organelle'][0] == 'mitochondrion':
            mol_type = 'm'
    if mol_type is None:
        if loci.get('keyless') and loci.get('keyless').get('Protein'):
            mol_type = 'p'

    return mol_type


def _extract_annotations(record):
    """
    Extract record annotations.

    :arg SeqRecord record: A biopython record.
    :return: Annotations as a dictionary.
    :rtype: dict
    """
    annotations = {}
    for annotation in record.annotations:
        if annotation in ANNOTATIONS:
            annotations[annotation] = record.annotations[annotation]

    return annotations


def _build_info(record, reference, loci):

    info = _extract_annotations(record)

    if reference.source.orientation:
        info['orientation'] = reference.source.orientation
    info.update(reference.source.qualifiers)

    info['mol_type'] = _extract_mol_type(reference, loci)

    info['id'] = '{}.{}'.format(
        info.get('accessions')[0],
        info.get('sequence_version'))

    return info


def _add_extra_non_feature_loci(loci, protein_id):
    """
    For NPs the mRNA and gene features need to be artificially constructed.
    """
    gene_keys = ['gene', 'gene_synonym', 'HGNC', 'GeneID', 'CCDS']
    gene = {}
    mrna = {}
    cds = {}

    if loci.get('keyless') and loci.get('keyless').get('CDS'):
        if len(loci.get('keyless').get('CDS')) == 1:
            cds = loci.get('keyless').get('CDS')[0]
            cds.qualifiers['protein_id'] = protein_id
            gene_name = cds.qualifiers.get('gene')
            if gene_name and (gene_name not in gene):
                gene_locus = Locus()
                gene_locus.qualifiers = {k: v for k, v in cds.qualifiers.items()
                                         if k in gene_keys}
                gene[gene_name] = gene_locus
            mrna_id = cds.qualifiers.get('coded_by').split(':')[0]
            if mrna_id and (mrna_id not in mrna):
                mrna_locus = Locus()
                mrna_locus.locus_type = 'mRNA'
                mrna_locus.qualifiers = {k: v for k, v in cds.qualifiers.items()
                                         if k in gene_keys}
                mrna_locus.qualifiers['transcript_id'] = mrna_id
                mrna_locus.qualifiers['protein_id'] = protein_id
                mrna[mrna_id] = mrna_locus
        else:
            return

    if gene:
        if loci.get('gene'):
            loci['gene'].update(gene)
        else:
            loci['gene'] = gene

    if mrna:
        if loci.get('mRNA'):
            loci['mRNA'].update(mrna)
        else:
            loci['mRNA'] = mrna
    if cds:
        if loci.get('CDS'):
            loci['CDS'].update({protein_id: cds})
        else:
            loci['CDS'] = {protein_id: cds}


def _construct_dependencies(loci):
    """
    Add loci to their parents (e.g., mRNAs to genes) and link them to their
    pairs (e.g., mRNAs to CDSs).

    :arg dict loci: Loci dictionary.
    """
    if loci.get('gene') and loci.get('mRNA') and loci.get('CDS'):
        for mrna in loci['mRNA']:
            cds_link = link_reference(mrna)[0]
            if cds_link and (cds_link in loci['CDS']):
                loci['CDS'][cds_link].link = loci['mRNA'][mrna]
                loci['mRNA'][mrna].link = loci['CDS'][cds_link]
                print('added', mrna)
            mrna_gene = loci['mRNA'][mrna].qualifiers.get('gene')
            if mrna_gene and mrna_gene in loci['gene']:
                loci['gene'][mrna_gene].add_child(loci['mRNA'][mrna])
        for cds in loci['CDS']:
            cds_gene = loci['CDS'][cds].qualifiers.get('gene')
            if cds_gene and cds_gene in loci['gene']:
                loci['gene'][cds_gene].add_child(loci['CDS'][cds])


def print_loci(loci):
    """
    Simple loci print for debug purposes.

    :param loci: Loci reference model instance.
    """
    for gene in loci['gene']:
        print('{}:'.format(loci['gene'][gene].qualifiers.get('gene')))
        if 'mRNA' in loci['gene'][gene].children:
            for child in loci['gene'][gene].children['mRNA']:
                print(' {} - {}'.format(
                    child.qualifiers.get('transcript_id'),
                    child.link.qualifiers.get('protein_id')))


