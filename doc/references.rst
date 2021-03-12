Transcript references
#####################

Features that should be included:

RefSeq ones, which are curated, contain only one transcript per reference.

In a corresponding genbank file there should be the following features:

- ``source``: only one.
- ``gene``: only one.
- ``CDS``: only one.
- ``exon``: multiple.

Important things:

- The ``protein_id`` and ``product`` appear in the ``CDS`` feature.
- The exon locations should be constructed from the ``exon`` features.
- For RefSeq references the ``transcript_id`` is the reference file name and
  should be gathered from there. Alternatively,

To be checked:

- If the exons are one after the other.

To be derived/constructed:

- The ``transcript_id``: from the reference name (alternatively, from the
  annotation section).
- The ``transcript_name``, from the ``gene`` name in the qualifiers.

Questions:

- Should we infer the CDS features?

Strange transcript references
=============================

AA010203.1
^^^^^^^^^^

Contains only the ``source`` feature.

NM_152263.2
^^^^^^^^^^^

Contains no ``exon`` features. We should mention this in an error message.

L41870.1
^^^^^^^^

Protein references
##################


Genomic references
##################
