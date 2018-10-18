"""
Main gbparser module.

@requires: biopython
"""
import os
import logging
import io
import json
import hashlib

from Bio import SeqRecord
from Bio import SeqFeature
from Bio import SeqIO


class PositionError(Exception):
    """
    Exception raised for errors in the position input.
    """
    pass


class Position:
    """
    A position in the sequence.

    At the moment the following positions are supported:
    - exact positions: a position which is specified as exact along the
      sequence and is represented as a number.
    - before positions: a fuzzy position that occurs prior to some coordinate.
      Example:
          `<20' - the real position is located somewhere less than 20.
    - after positions: a fuzzy position that occurs after to some coordinate.
      Example:
          `>20' - the real position is located somewhere after 20.

    The maximum position is determined by the int size, which could suffice
    for an entire chromosome (Human Chromosome 1 is around 250M bases).

    Assuming p1, p1 as Position instances and c as int, the following
    operations are allowed:
        - p1 + c
        - c + p1
        - p1 - c
        - c - p1
        - p1 - p2
    Note that p1 + p2 is not allowed. There is no reason at the moment to add
    two positions.
    """
    # TODO: Decide exactly what operations to allow for positions.
    # TODO: Allow only positions > 0.
    def __init__(self, position=None):
        """
        The position initializer.

        :param position: A position location as int, str, or Position.
        """
        self._position = None
        self._fuzzy = None

        if position is not None:
            self.position = position

    @property
    def position(self):
        """
        Getter for the position value.
        Note that the position value is returned for both fuzzy and exact
        positions. Check position type with is_fuzzy() method.

        :return: Position value.
        """
        # if self.is_fuzzy():
        #     warnings.warn('Not exact position.')
        return self._position

    @position.setter
    def position(self, position):
        """
        Position setter.

        :param position: A position location as int, str, or Position.
        """
        # If supplied position is an int then it is an exact position.
        if isinstance(position, int):
            self._position_from_int(position)
        # If supplied position is a string, then we should
        # determine if it is an exact or a fuzzy position.
        elif isinstance(position, str):
            self._position_from_str(position)
        # If supplied position is another Position instance we
        # assume that the string representation is correct.
        elif isinstance(position, Position):
            self._position_from_str(str(position))
        # If position is an instance of something else we raise an error.
        else:
            raise PositionError('Incorrect position input.')

    @property
    def fuzzy(self):
        """
        Getter for the fuzzy description.

        :return: Position value.
        """
        return self._fuzzy

    @fuzzy.setter
    def fuzzy(self, fuzzy):
        """
        Setter for the fuzzy positions.

        :param fuzzy: Fuzzy position type '<' or '>'.
        """
        if fuzzy in ['<', '>']:
            self._fuzzy = fuzzy

    def _position_from_int(self, position):
        """
        An exact position is created by providing an integer value.
        In this case the fuzzy position attribute is set to None.

        :param position: Position value as int.
        """
        self._position = position
        self._fuzzy = None

    def _position_from_str(self, position):
        """
        By providing a string input both fuzzy (before or after) and exact
        positions can be created.

        :param position: Position value as string.
        :return:
        """
        # Check if position is empty.
        if not position:
            raise PositionError('Incorrect position input.')
        # Check for fuzzy before/after positions.
        if position[0] in ['<', '>']:
            try:
                pos = int(position[1:])
                self._position = pos
                self._fuzzy = position[0]
            except ValueError:
                raise PositionError('Incorrect position input.')
            return
        # Check for exact position.
        # With a string you can also define an exact position
        # if '<' or '>' are missing.
        try:
            pos = int(position)
            self._position_from_int(pos)
        except ValueError:
            raise PositionError('Incorrect position input.')

    def is_fuzzy(self):
        """
        Checks if a position is fuzzy or not.

        :return: True of position is fuzzy, False otherwise.
        """
        if self._fuzzy is None:
            return False
        else:
            return True

    def is_before(self):
        """
        Checks if a position is of 'before' type, e.g., '<100'.

        :return: True of position is a 'before' one or False otherwise.
        """
        return bool(self._fuzzy == '<')

    def is_after(self):
        """
        Checks if a position is of 'after' type, e.g., '>100'.

        :return: True of position is an 'after' one or False otherwise.
        """
        return bool(self._fuzzy == '>')

    def add_int(self, value):
        """
        Adds the specified amount to the position.

        :param value: Amount to be added.
        """
        if isinstance(value, int):
            self._position += value
        else:
            raise PositionError('Different argument than int provided.')

    def __str__(self):
        output = ''
        if self.is_fuzzy():
            output += self._fuzzy
        output += str(self._position)
        return output

    def __add__(self, other):
        if isinstance(other, int):
            return self._position + other
        return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            return self._position - other
        if isinstance(other, Position):
            return self._position - other
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, int):
            return other - self._position


class Locus(object):
    """
    A Locus object, to store data about the genes, mRNA, CDS, etc. features.
    """
    def __init__(self, start=None, end=None, orientation=None,
                 locus_type=None, parts=None, sequence=None,
                 qualifiers=None, parent=None, children=None,
                 config=None, link=None):
        """
        :param qualifiers: Dictionary containing the genbank qualifiers.
                           Taken directly from the BioPython SeqFeature.
        :param start: Locus beginning. A position location.
        :param end: Locus end. A position location.
        :param type: Locus type, e.g., exon, CDS, mRNA, etc..
        :param parts: The composing parts for compound sequences.
        :param orientation: Locus orientation:
                            - *+1* for *3'* to *5'* strand orientation
                            - *-1* for *5'* to *3'* strand orientation.
        :param sequence: The actual nucleotide bases, protein, etc. sequence.
        """
        self._start = Position(start)
        self._end = Position(end)
        # Check that *start* is greater than the *end*.
        if None not in [start, end]:
            if self._start.position > self._end.position:
                raise Exception('Locus start > sequence end.')

        self._orientation = None
        if orientation is not None:
            self.orientation = orientation

        self._locus_type = None
        if locus_type is not None:
            self.locus_type = locus_type

        self._parts = []
        if parts is not None:
            self.parts = parts

        self._sequence = None
        if sequence is not None:
            self.sequence = sequence

        self._qualifiers = {}
        if qualifiers is not None:
            self.qualifiers = qualifiers

        self._parent = None
        if parent is not None:
            self.parent = parent

        self._children = {}
        if children is not None:
            self.children = children

        self._config = {}
        if config is not None:
            self.config = config

        self._link = {}
        if link is not None:
            self.link = link

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, position):
        self._start = position

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, position):
        self._end = position

    @property
    def locus_type(self):
        return str(self._locus_type)

    @locus_type.setter
    def locus_type(self, seq_type):
        self._locus_type = seq_type

    @property
    def orientation(self):
        return self._orientation

    @orientation.setter
    def orientation(self, orientation):
        # Check if orientation is correct.
        if orientation not in [None, 1, -1]:
            raise Exception('Incorrect sequence orientation (not 1 or -1).')
        self._orientation = orientation

    @property
    def parts(self):
        return self._parts

    @parts.setter
    def parts(self, parts):
        # *parts* should be a list of Locus instances, so we check that.
        if isinstance(parts, list):
            # Check if providing parts are all Loci.
            for part in parts:
                if not isinstance(part, Locus):
                    raise Exception('Incorrect sequence part element provided.')
            self._parts = parts
        # But *parts* can be also just a Locus instance. If this is the
        # case we just create the list with that instance in it.
        elif isinstance(parts, Locus):
            self.parts = [].append(parts)
        else:
            raise Exception('Incorrect parts provided.')

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence

    @property
    def qualifiers(self):
        return self._qualifiers

    @qualifiers.setter
    def qualifiers(self, qualifiers):
        self._qualifiers = qualifiers

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        self._parent = parent

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, children):
        self._children = children

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, config):
        self._config = config

    @property
    def link(self):
        return self._link

    @link.setter
    def link(self, link):
        self._link = link

    def __len__(self):
        """
        Get the sequence length. We assume that positions are positive and
        that start < end. Orientation could be implemented in future.

        :return: The sequence length.
        """
        # TODO: Check if negative positions or orientation is provided.
        if None not in [self._end.position, self._start.position]:
            return self._end - self._start + 1
        else:
            raise Exception('None limits present in sequence.')

    def __str__(self):
        """
        Simple sequence information string representation.
        """
        output = ''
        # if self.get_qualifier(self.config['key']):
        #     output += 'sequence key: %s\n' % self.qualifiers[self.config['key']]
        output += 'start: %s\n' % str(self._start)
        output += 'end: %s\n' % str(self._end)
        output += json.dumps(self.get_parts_string())
        output += '\nstrand: %s\n' % self.get_strand()
        return output

    def to_dict(self):
        """
        Converts the sequence into a dictionary. It includes the children as
        well. Called in principle to help jsonify the record.

        :return: The sequence as a dictionary.
        """
        values = dict(start=str(self.start),
                      end=str(self.end),
                      strand=self.get_strand(),
                      parts=[])
        for part in self.parts:
            values['parts'].append(part.to_dict())

        values.update(self.qualifiers)

        if 'children' in self.config:
            for child_type in self.children:
                values[child_type] = []
                for child in self.children[child_type]:
                    values[child_type].append(child.to_dict())

        if 'key' in self.config:
            if self.get_qualifier(self.config['key']):
                output = {self.qualifiers[self.config['key']]: values}
            else:
                output = {'None': values}
        else:
            output = values
        return output

    def fuzzy_position_inside(self):
        """
        To detect whether a fuzzy position is present in the start/stop
        positions of the sequence, as well as in any start/stop positions of
        its composing parts.
        :return: True if fuzzy position found in sequence, False otherwise.
        """
        if self.start.is_fuzzy():
            return True
        if self.end.is_fuzzy():
            return True
        for part in self.parts:
            if part.start.is_fuzzy():
                return True
            if part.end.is_fuzzy():
                return True
        return False

    def get_qualifier(self, key):
        """
        Checks if the provided 'key' is present in the '_qualifiers'
        dictionary and returns the value if so, or None otherwise.

        :param key: Key to be searched for.
        :return: The value for the specified key if found, or None.
        """
        if self._qualifiers is not None:
            if isinstance(self._qualifiers, dict):
                if str(key) in self._qualifiers:
                    return self._qualifiers[str(key)]
        return None

    def add_qualifier(self, key, value):
        """
        Adds the provided qualifier to '_qualifiers' dictionary.

        :param key: provided key.
        :param value: provided value.
        """
        self._qualifiers[key] = value

    def update_qualifier(self, key, value):
        """
        Update the provided qualifier to '_qualifiers' dictionary.

        :param key: provided key.
        :param value: provided value.
        """
        self._qualifiers[key] = value

    def in_gene(self):
        """
        Checks if the Locus is part of a gene or not. It returns the name
        of the gene if so and None otherwise.
        Note that for this it only checks for the appearance of 'gene' as key
        in the 'qualifiers' dictionary attribute. It can still be part of a
        gene if the 'gene' key is missing, but this should be checked in a
        different manner.

        :return: The name of the gene or 'None' if the gene was not found.
        """
        # Checking first if the '_qualifiers' attribute is not None.
        if self._qualifiers is not None:
            # '_qualifiers' should be also a dictionary.
            if isinstance(self._qualifiers, dict):
                if 'gene' in self._qualifiers:
                    # Gene found, so we return it.
                    # It seems like BioPyothon provides stores the gene
                    # name in a list, so we should get the first element.
                    return self._qualifiers['gene'][0]
        # No 'gene' found among the '_qualifiers' dictionary.
        return None

    def get_strand(self):
        """
        Converts the orientation into either '+' or '-'.
        """
        if self._orientation == -1:
            return '-'
        if self._orientation == 1:
            return '+'

    def add_child(self, seq):
        """
        Adds a child to the children list.
        :param seq: child sequence to be added.
        """
        if seq.locus_type in self.children:
            self.children[seq.locus_type].append(seq)
        else:
            self.children[seq.locus_type] = [seq]

    def get_child(self, feature_type, q_k, q_v):
        """
        Accessor for a specific child in the children list.

        :param feature_type: the feature of the child.
        :param q_k: the qualifier key to look for.
        :param q_v: the qualifier value to look for.
        """
        children_of_type = self.children.get(feature_type)
        if children_of_type:
            for child in children_of_type:
                if child.get_qualifier(q_k) == q_v:
                    return child

    def append_part(self, part):
        """
        Append a composing Locus to the '_parts' attribute.

        :param part: A Locus instance.
        """
        # Create the '_parts' list if None.
        if self._parts is None:
            self._parts = []
        # Check if the provided part is a Locus instance.
        if isinstance(part, Locus):
            self._parts.append(part)
        else:
            raise Exception('Part to be appended is not a sequence.')

    def get_parts_string(self):
        """
        Returns the start and end positions as strings part of a dictionary.

        :return: Dictionary with start and end positions as string values.
        """
        starts = ''
        ends = ''
        delimiter = ''
        for part in self._parts:
            starts += delimiter + str(part.start)
            ends += delimiter + str(part.end)
            delimiter = ','
        return {'start_positions': starts, 'end_positions': ends}

    def get_key_type_and_value(self):
        """

        :return:
        """
        key_type = self.config.get('key')
        return key_type, self.qualifiers.get(key_type)


class Reference:
    """
    Structured container of all the features in a genbank file.
    """
    def __init__(self, start=None, end=None,
                 reference=None, cfg=None):
        # All the features present in the input file.
        self.features = {}

        self.loci = []

        # All the parents in a tree with feature types as keys.
        self._parents_dict = {}

        # All the sequences unattached to the record tree.
        self._unattached = []

        # All the sequence without the key present in the qualifiers.
        self._keyless_sequences = []

        # The number of sequences with fuzzy positions.
        self._fuzzy_sequences = 0

        # Reference content (usually the raw genbank file content).
        self._input_record = None

        # Annotations present in the record.
        self._annotations = {}

        if reference and cfg:
            self.create_record(reference, cfg)

    def child_of_type(self, feature):
        """
        Finds the parent type of the provided feature.
        :param feature: a configuration feature dictionary.
        :return:
        """
        if 'parent' in feature:
            if isinstance(feature['parent'], dict):
                if 'type' in feature['parent']:
                    return feature['parent']['type']
            else:
                return feature['parent']

    def child_of_key(self, feature):
        """
        Get the parent key with which the child should be added to the parent.
        :return: The parent key or None.
        """
        if 'parent' in feature:
            if isinstance(feature['parent'], dict):
                if 'key' in feature['parent']:
                    return feature['parent']['key']

    def create_parents_dict(self):
        """
        Initializes the parents dictionary with the keys from the configuration.
        :param features: dictionary with input fields according to the
                         configuration format.
        """
        if 'features' in self.config:
            features = self.config['features']
            for feature in features:
                parent_type = self.child_of_type(features[feature])
                parent_key = self.child_of_key(features[feature])
                if parent_type:
                    self._parents_dict[parent_type] = {} if parent_key else []
        # The record is added also among the parents.
        # This should be always performed.
        self._parents_dict['record'] = [self]

    def add_to_parent_with_key(self, child):
        """
        Tries to add a child sequence to a parent.

        :param child: The sequence to be added to the parent.
        :return True if addition performed and False otherwise.
        """
        parent_type = self.child_of_type(child.config)
        if parent_type and parent_type not in self._parents_dict:
            return
        parent_key = child.get_qualifier(self.child_of_key(child.config))

        if self.child_of_type(child.config) in self._parents_dict:
            parent_type = self.child_of_type(child.config)
            if parent_type and parent_type in self._parents_dict:
                if parent_key:
                    parent = self._parents_dict[parent_type][parent_key]
                    parent.add_child(child)
                    child.parent = parent
                    return True
        return False

    def add_to_parents_dict(self, seq):
        """

        :param seq:
        """
        parent_node = self._parents_dict.get(seq.locus_type)
        if isinstance(parent_node, dict):
            parent_type = self.child_of_key(seq.config)
            key_type = seq.config['key']
            key = seq.get_qualifier(key_type)
            parent_node[key] = seq
        if isinstance(parent_node, list):
            parent_node.append(seq)

    def extract_features(self, input_record):
        """
        Loops over the biopython record features and extracts the ones
        specified in the configuration file. The extracted feature is
        converted into a Locus instance.

        The Loci tree is formed at the same time, with parents being
        added to the parents_dict, if found, and children linked

        :param input_record: A biopython record.
        """
        for feature in self.config['features']:
            self._features[feature] = []
        for feature in input_record.features:
            # Check if the feature is in the configuration features.
            if feature.type in self.config['features']:
                seq = Locus()
                # Extracting the configuration for the sequence.
                configuration = self.config['features'][feature.type]
                # Transform the biopython feature into a Locus instance.
                seq.convert_biopython_feature(feature, configuration)
                # Check is the sequence contains fuzzy positions.
                if seq.fuzzy_position_inside():
                    self._fuzzy_sequences += 1
                # Finding if a sequence has a key with a value.
                key, value = seq.get_key_type_and_value()
                # We disregard the keyless sequences but we keep track of
                # their number, which is added to the log file later.
                if key and not value:
                    self._keyless_sequences.append(seq)
                else:
                    self._features[feature.type].append(seq)
                    # Try to add it to the _parents_dict.
                    self.add_to_parents_dict(seq)
                    # If sequence is a child then we add it to the parent.
                    parent = self.child_of_type(self.config['features'][feature.type])
                    if parent:
                        # We do not forget to add the record children.
                        if parent == 'record':
                            self.add_child(seq)
                        else:
                            # It seems that either the parent was not added yet
                            # (best case) or there is no parent.
                            if not self.add_to_parent_with_key(seq):
                                self._unattached.append(seq)

    def extract_annotations(self, input_record):
        """
        Extract the annotations fields from the record, according to the
        provided configuration.

        :param input_record: A biopython record.
        """
        for annotation in input_record.annotations:
            # Check if the feature is in the configuration features.
            if annotation in self.config['annotations']:
                self._annotations[annotation] = input_record.annotations[annotation]

    def attach_the_unattached(self):
        """
        Attaches all the remaining unattached sequence to their parents.
        Used primarily after running the extract_features method.
        """
        still_unattached = []
        for seq in self._unattached:
            if self.add_to_parent_with_key(seq):
                still_unattached.append(seq)
        self._unattached = still_unattached
        if len(still_unattached) > 0:
            logging.info('Still %s children unattached.', len(still_unattached))
        else:
            logging.info('All children attached.')

    def perform_linking(self):
        """
        Call any linking functions that appear in the configuration.
        """
        if self.config.get('link') is None:
            return
        logging.info('Link started.')
        for link in self.config['link']:
            link_main_function = 'link_' + link
            if self.__getattribute__(link_main_function):
                self.__getattribute__(link_main_function)(self.config['link'][link])
            else:
                logging.info('No %s function implemented.', link_main_function)
        logging.info('Link ended.')

    def link_mrna_with_cds(self, cfg):
        """
        Calls the corresponding function to perform the actual linking of an
        mrna sequence to a cds sequence. Such a function should exist. All
        these should be actually placed in a separate file.

        :param cfg: the configuration (should contain the link function type
                    should be present).
        """
        for link_function_type in cfg:
            link_function = 'link_' + link_function_type
            if self.__getattribute__(link_function):
                self.__getattribute__(link_function)(cfg[link_function_type])

    def link_by_ncbi_file_website(self, cfg):
        """
        Link the gene transcripts with proteins via the ncbi website.
        The corresponding .gb file of the transcript is downloaded and the
        protein_id is searched in the CDS feature of the record file.

        :param cfg:
        """
        for gene in self._parents_dict['gene']:
            mrnas = self._parents_dict['gene'][gene].children.get('mRNA')
            cdss = self._parents_dict['gene'][gene].children.get('CDS')
            if mrnas:
                for mrna in mrnas:
                    transcript_id = mrna.get_qualifier('transcript_id')
                    # logging.info('Linking transcript id %s', str(transcript_id))
                    # try:
                    #     p_accession = transcript_to_protein2(
                    #         transcript_id, cfg['options']['save_file'])
                    # except Exception as exception:
                    #     logging.info('Error during transcript to protein link:'
                    #                  ' "%s" Transcript id: "%s"',
                    #                  str(exception), str(transcript_id))
                    # else:
                    #     mrna.update_qualifier('protein_id_link', p_accession)
                    #     if isinstance(cdss, list):
                    #         for cds in cdss:
                    #             if p_accession == cds.get_qualifier('protein_id'):
                    #                 mrna.link = cds
                    #                 cds.link = mrna
                    #                 cds.update_qualifier('transcript_id_link',
                    #                                      transcript_id)
                    #                 mrna.update_qualifier('cds_start',
                    #                                       str(cds.start))
                    #                 mrna.update_qualifier('cds_end',
                    #                                       str(cds.end))
                    #                 mrna.update_qualifier('locus_tag',
                    #                                       cds.qualifiers.get('locus_tag'))
                    #                 mrna.update_qualifier('codon_start',
                    #                                       cds.qualifiers.get('codon_start'))
                    #                 # mrna.update_qualifier('cds_parts', cds.parts)

    def link_by_mapping_file(self, cfg):
        """
        Link gene transcripts with proteins via a provided mapping file.
        The content of the file should utilize the following format:

            accession.version, accession.version

        with the first "accession.version" begin for the transcript, and the
        second one for the protein.

        Configuration dictionary fields:
        - path: the path of the file

        :param cfg: Provided function configuration.
        :return:
        """
        logging.info('Linking by mapping file.')
        # Load the file
        mapping_file = cfg['path']
        if os.path.isfile(mapping_file):
            logging.info('Path to mapping file: %s', mapping_file)
            mappings = {}
            with open(mapping_file, 'r') as mapping_file:
                for line in mapping_file.readlines():
                    transcript_id = line.strip().split(',')[0]
                    protein_id = line.strip().split(',')[1].strip()
                    mappings[transcript_id] = protein_id
        else:
            logging.info('No mapping file found.')
            return
        unlinked = 0
        for gene in self._parents_dict['gene']:
            mrnas = self._parents_dict['gene'][gene].children.get('mRNA')
            cdss = self._parents_dict['gene'][gene].children.get('CDS')
            if mrnas:
                for mrna in mrnas:
                    linked = False
                    transcript_id = mrna.get_qualifier('transcript_id')
                    protein_id = mappings.get(str(transcript_id))
                    if protein_id:
                        mrna.update_qualifier('protein_id_link', protein_id)
                        if isinstance(cdss, list):
                            for cds in cdss:
                                if protein_id == cds.get_qualifier('protein_id'):
                                    mrna.link = cds
                                    cds.link = mrna
                                    cds.update_qualifier('transcript_id_link',
                                                         transcript_id)
                                    mrna.update_qualifier('cds_start',
                                                          str(cds.start))
                                    mrna.update_qualifier('locus_tag',
                                                          cds.qualifiers.get('locus_tag'))
                                    mrna.update_qualifier('codon_start',
                                                          cds.qualifiers.get('codon_start'))
                                    # mrna.update_qualifier('cds_parts', cds.parts)
                                    linked = True
                    else:
                        logging.info('No protein_id found in mapping file for'
                                     'transcript_id %s', transcript_id)
                    if not linked:
                        # print('\nunlinked')
                        unlinked += 1
        logging.info('Unlinked transcripts: %s', unlinked)

    def get_biopython_record(self, reference):
        """
        Transforms the reference into a biopython record and adds the checksum
        to the qualifiers dictionary.

        :param reference: accession.version or file path.
        """
        if os.path.isfile(reference):
            logging.info('The reference is a file.')
            hash_md5 = hashlib.md5()
            content = ''
            with open(reference, 'r') as file:
                for line in file:
                    hash_md5.update(line.strip().encode('utf-8'))
                    content += line
            checksum_reference = hash_md5.hexdigest()
        else:
            logging.info('Reference is in the format accession.version - retrieve from NCBI.')
            # try:
            #     content, checksum_reference = get_reference(reference)
            # except Exception as exception:
            #     logging.info('Reference %s not retrieved due to error %s',
            #                  reference, str(exception))
        self._input_record = content
        self.qualifiers['checksum_reference'] = checksum_reference
        input_record = io.StringIO(content)
        try:
            biopython_record = SeqIO.read(input_record, 'genbank')
        except:
            logging.info('Content for reference %s could not be read by biopython.', reference)
            return None
        else:
            return biopython_record

    def write_json(self, dir_path):
        """
        Writes the json formatted record inside the provided directory path.
        The file name is 'accession_version.json'.

        :param dir_path: Directory in which the record is saved.
        """
        output_path = dir_path + '/' + self.get_accession_version() + '.json'
        try:
            record_json = json.loads(self.to_json())
        except Exception as exception:
            logging.debug('Record load to json failed with %s', str(exception))
            return
        try:
            with open(output_path, 'w') as output_file:
                output_file.write(json.dumps(record_json, indent=2))
        except Exception as exception:
            logging.debug('Reference %s Json file writing error: %s',
                          self.reference, str(exception))
        else:
            logging.info('Json saved to: %s', output_path)

    def write_input(self, dir_path):
        """
        Writes the record input to a file
        :param dir_path:
        """
        output_path = dir_path + '/' + self.get_accession_version() + '.gb'
        try:
            with open(output_path, 'w') as output_file:
                output_file.write(self._input_record)
        except Exception as exception:
            logging.debug('Write input to file failed with: %s', str(exception))
        else:
            logging.info('Input saved to: %s', output_path)

    def write_sequence(self, dir_path):
        checksum_sequence = self.qualifiers.get('checksum_sequence')
        if checksum_sequence is None:
            logging.info('Locus without checksum - not saved.')
            return

        output_path = dir_path + '/' + checksum_sequence + '.sequence'
        try:
            output_file = open(output_path, 'w')
        except IOError as exception:
            logging.debug('Locus write to file failed with %s', str(exception))
        else:
            output_file.write(self.sequence)
            logging.info('Locus saved to: %s', output_path)

    def save_output(self):
        """
        If required (mentioned in the configuration), it saves the output to
        the files specified in the configuration. This may include the record
        in json format, the input genbank file (from which the record was
        extracted), and the actual sequence.
        """
        if self.config.get('output') is None:
            return

        cfg_output = self.config['output']
        if not isinstance(cfg_output, dict):
            logging.info('Provided output configuration not a dictionary.')
            return
        dir_path = cfg_output.get('path')
        if not dir_path:
            logging.info('No output path in configuration.')
            return
        save_list = cfg_output.get('save')
        if not save_list or not isinstance(save_list, list):
            return
        if 'json' in save_list:
            self.write_json(dir_path)
        if 'input' in save_list:
            self.write_input(dir_path)
        if 'sequence' in save_list:
            self.write_sequence(dir_path)

    def create_record(self, reference, cfg):
        """
        Creates a gbparser record from a genbank reference file according to
        the provided configuration.

        The biopython record is first obtained.

        :param reference: The reference (file path or only 'accession[.version]').
        :param cfg: The configuration based on which the record is created.
        """
        input_record = self.get_biopython_record(reference)
        # A provided input_record must be an instance of BioPython SeqRecord.
        if not isinstance(input_record, SeqRecord.SeqRecord):
            logging.info('Record_input not of BioPython SeqRecord type.')
            raise Exception('Record_input not of BioPython SeqRecord type.')
        # The provided input_record must contain some features, otherwise
        # we can stop.
        if not input_record.features:
            return

        self.config = cfg
        self.create_parents_dict()
        if self.config['features']:
            self.extract_features(input_record)
        if self.config['annotations']:
            self.extract_annotations(input_record)
        self.attach_the_unattached()
        self.perform_linking()
        # Extract the sequence, compute its checksum, and add it to qualifiers.
        self.sequence = str(input_record.seq)
        checksum_sequence = hashlib.md5(self.sequence.encode('utf-8')).hexdigest()
        self.qualifiers['checksum_sequence'] = checksum_sequence
        self.qualifiers.update(self.get_from_source())
        self.save_output()
        self.qualifiers['transl_table'] = ''
        self.qualifiers['source'] = 'NCBI eutils'

    def link_transcript_with_protein(self):
        """
        Tries to link all the record transcripts with the proteins.
        """
        for mrna in self._sequences['mRNA']:
            transcript_id = mrna.get_qualifier('transcript_id')[0]
            accession = str(transcript_id).split('.')[0]
            print(transcript_id, end='\r')
            # try:
            #     p_accession = transcript_to_protein2(transcript_id)
            # except Exception as exception:
            #     print('Error during transcript to protein link:\n ',
            #           str(exception), 'Transcript id: ' + transcript_id)
            # else:
            #     mrna.update_qualifier('protein_id_link', p_accession)

    def get_accession(self):
        """
        Retrieves the record accession present in the annotations.

        :return: record accession as string.
        """
        if self._annotations and 'accessions' in self._annotations:
            if isinstance(self._annotations['accessions'], list) and \
                            len(self._annotations) > 0:
                return str(self._annotations['accessions'][0])

    def get_version(self):
        """
        Retrieves the record version number present in the annotations.

        :return: record version number as string.
        """
        if self._annotations and 'sequence_version' in self._annotations:
            return str(self._annotations['sequence_version'])

    def get_accession_version(self):
        """
        Retrieves both record accession and version number from the annotations.
        :return: record 'accession.version' as string.
        """
        accession = self.get_accession()
        version = self.get_version()
        if accession and version:
            return str(accession) + '.' + str(version)

    def get_transcripts(self):
        """
        Yields all the transcripts in order to be added to the DB.

        :return:
        """
        pass

    def get_details(self):
        source = self.children.get('source')[0]
        if source:
            details = {'accession': self.get_accession(),
                       'version': self.get_version(),
                       'mol_type': source.get_qualifier('mol_type'),
                       'start': str(source.start),
                       'stop': str(source.end),
                       'organism': source.get_qualifier('organism')}
            return details

    def get_from_source(self):
        source = self.children.get('source')[0]
        if source:
            return {
                'accession': self.get_accession(),
                'version': self.get_version(),
                'date_annotation': self._annotations.get('date'),
                'length': int(source.end - source.start + 1),
                'mol_type': source.get_qualifier('mol_type'),
            }

    def get_from_cds(self):
        cds = self.children.get('gene').get('cds')[0]
        if cds:
            return {
                'transl_table': cds.get_qualifier('transl_table'),
            }

    def to_json(self):
        """
        Converts the record into JSON.

        :return: record JSON representation (string).
        """
        output = {'annotations': self._annotations,
                  'features': {}}
        for child_type in self.children:
            output['features'][child_type] = []
            for child in self.children[child_type]:
                if child_type == child.locus_type:
                    output['features'][child_type].append(child.to_dict())
        return json.dumps(output, indent=2)

    def __str__(self):
        """
        Simple record string representation.
        """
        return 'Total loci: {}'.format(len(self.loci))
        # return json.dumps(self._annotations)
