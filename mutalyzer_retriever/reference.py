"""
Main gbparser module.

@requires: biopython
"""
import hashlib
import io
import json
import logging
import os

from Bio import SeqFeature, SeqIO, SeqRecord


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
            raise PositionError("Incorrect position input.")

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
        if fuzzy in ["<", ">"]:
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
            raise PositionError("Incorrect position input.")
        # Check for fuzzy before/after positions.
        if position[0] in ["<", ">"]:
            try:
                pos = int(position[1:])
                self._position = pos
                self._fuzzy = position[0]
            except ValueError:
                raise PositionError("Incorrect position input.")
            return
        # Check for exact position.
        # With a string you can also define an exact position
        # if '<' or '>' are missing.
        try:
            pos = int(position)
            self._position_from_int(pos)
        except ValueError:
            raise PositionError("Incorrect position input.")

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
        return bool(self._fuzzy == "<")

    def is_after(self):
        """
        Checks if a position is of 'after' type, e.g., '>100'.

        :return: True of position is an 'after' one or False otherwise.
        """
        return bool(self._fuzzy == ">")

    def add_int(self, value):
        """
        Adds the specified amount to the position.

        :param value: Amount to be added.
        """
        if isinstance(value, int):
            self._position += value
        else:
            raise PositionError("Different argument than int provided.")

    def __str__(self):
        output = ""
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

    def __init__(
        self,
        start=None,
        end=None,
        orientation=None,
        locus_type=None,
        parts=None,
        sequence=None,
        qualifiers=None,
        parent=None,
        children=None,
        config=None,
        link=None,
    ):
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
                raise Exception("Locus start > sequence end.")

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
            raise Exception("Incorrect sequence orientation (not 1 or -1).")
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
                    raise Exception("Incorrect sequence part element provided.")
            self._parts = parts
        # But *parts* can be also just a Locus instance. If this is the
        # case we just create the list with that instance in it.
        elif isinstance(parts, Locus):
            self.parts = [].append(parts)
        else:
            raise Exception("Incorrect parts provided.")

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
        Get the locus length. We assume that positions are positive and
        that start < end. Orientation could be implemented in future.

        :return: The locus sequence length.
        """
        # TODO: Check if negative positions or orientation is provided.
        if None not in [self._end.position, self._start.position]:
            return self._end - self._start + 1
        else:
            raise Exception("None limits present in locus.")

    def __str__(self):
        """
        Simple locus information string representation.
        """
        if self._locus_type:
            feature_location = "{:10}{}..{}\n".format(
                self._locus_type, self._start, self._end
            )
        else:
            feature_location = "  {}..{}\n".format(self._start, self._end)
        qualifiers = "\n".join(
            "  {:15}: {}".format(k, v) for k, v in self.qualifiers.items()
        )
        parts = self.get_parts_string()
        output = "{} qualifiers:\n{}\n strand: {}".format(
            feature_location, qualifiers, self.get_strand()
        )
        if parts != "":
            output += "\n parts:{}".format(parts)
        if self.children:
            output += " children:\n"
            for child in self.children:
                output += "{}".format(child.locus_type)
        return output

    def to_dict(self):
        """
        Converts the locus into a dictionary. It includes the children as
        well. Called in principle to help jsonify the record.

        :return: The locus as a dictionary.
        """
        values = dict(
            start=str(self.start), end=str(self.end), strand=self.get_strand(), parts=[]
        )
        for part in self.parts:
            values["parts"].append(part.to_dict())

        values.update(self.qualifiers)

        if "children" in self.config:
            for child_type in self.children:
                values[child_type] = []
                for child in self.children[child_type]:
                    values[child_type].append(child.to_dict())

        if "key" in self.config:
            if self.get_qualifier(self.config["key"]):
                output = {self.qualifiers[self.config["key"]]: values}
            else:
                output = {"None": values}
        else:
            output = values
        return output

    def fuzzy_position_inside(self):
        """
        To detect whether a fuzzy position is present in the start/stop
        positions of the locus, as well as in any start/stop positions of
        its composing parts.
        :return: True if fuzzy position found in locus, False otherwise.
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
                if "gene" in self._qualifiers:
                    # Gene found, so we return it.
                    # It seems like BioPyothon provides stores the gene
                    # name in a list, so we should get the first element.
                    return self._qualifiers["gene"][0]
        # No 'gene' found among the '_qualifiers' dictionary.
        return None

    def get_strand(self):
        """
        Converts the orientation into either '+' or '-'.
        """
        if self._orientation == -1:
            return "-"
        if self._orientation == 1:
            return "+"

    def add_child(self, locus):
        """
        Adds a child to the children list.
        :param locus: child locus to be added.
        """
        if locus.locus_type in self.children:
            self.children[locus.locus_type].append(locus)
        else:
            self.children[locus.locus_type] = [locus]

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
            raise Exception("Part to be appended is not a locus.")

    def get_parts_string(self):
        """
        Returns the start and end positions as strings part of a dictionary.

        :return: Dictionary with start and end positions as string values.
        """
        output = ""
        delimiter = "\n"
        for part in self._parts:
            output += "{}  {}..{}".format(delimiter, part.start, part.end)
        return output

    def get_parts_list(self):
        """
        Returns the start and end positions in a list.

        """
        output = []
        if self._parts:
            for part in self._parts:
                output.append(int(str(part.start)))
                output.append(int(str(part.end)))
        else:
            output = [str(self.start), str(self.end)]
        return output

    def get_key_type_and_value(self):
        """

        :return:
        """
        key_type = self.config.get("key")
        return key_type, self.qualifiers.get(key_type)


class Reference:
    """
    Structured container of all the features in a reference file.
    """

    def __init__(self, start=None, end=None, reference=None, cfg=None):

        self.type = None

        self.info = {}

        self.source = {}

        self.loci = {}

        self.mol_type = None

        # All the locis without the key present in the qualifiers.
        self.keyless_loci = []

        # All the parents in a tree with feature types as keys.
        self._parents_dict = {}

        # All the sequences unattached to the record tree.
        self._unattached = []

        # The number of sequences with fuzzy positions.
        self._fuzzy_sequences = 0

        # Reference content (usually the raw genbank file content).
        self._input_record = None

        # Annotations present in the record.
        self._annotations = {}

        if reference and cfg:
            self.create_record(reference, cfg)

    def to_dict(self):
        """
        Draft loci object serializer.

        :arg dict loci:
        """
        json_model = {"reference": self.info}
        loci_json = []

        if "curated" in self.loci:
            loci_json.append(self.loci["curated"])
        else:
            for gene_name, gene in self.loci["gene"].items():
                if "mRNA" in gene.children:
                    for child in gene.children["mRNA"]:
                        if not child.fuzzy_position_inside():
                            locus_json = {
                                "transcript_id": child.qualifiers.get("transcript_id"),
                                "HGNC": self.loci["gene"][gene_name].qualifiers.get(
                                    "HGNC"
                                ),
                                "gene": gene_name,
                            }
                            if child.orientation:
                                locus_json.update({"orientation": child.orientation})
                            if child.parts:
                                locus_json.update({"exons": child.get_parts_list()})
                            if child.link:
                                locus_json.update(
                                    {
                                        "location": [
                                            int(str(child.link.start)),
                                            int(str(child.link.end)),
                                        ],
                                        "protein_id": child.link.qualifiers.get(
                                            "protein_id"
                                        ),
                                    }
                                )

                            loci_json.append(locus_json)

        json_model.update({"loci": loci_json})
        return json_model

    def to_json(self):
        return json.dumps(self.to_dict(), sort_keys=True, indent=2)

    def __str__(self):
        """
        Simple record string representation.
        """
        # print(self.loci.values())
        info = "\n".join([" {:17}: {}".format(k, v) for k, v in self.info.items()])
        loci = "\n".join(map(str, self.loci.values()))

        return "\nReference info:\n{}\n{}\nTotal loci: {}".format(
            info, loci, len(self.loci)
        )
        # return json.dumps(self._annotations)
