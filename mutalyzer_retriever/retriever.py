from . import configuration, parser
from .sources import ensembl, lrg, ncbi


def fetch_annotations(reference_id, reference_type=None):
    """

    :arg str reference_id: The reference ID.
    :arg str reference_type: The source from where to retrieve the annotations.

    :returns tuple: Annotations, reference type, reference source.
    """
    annotations = lrg.fetch_lrg(reference_id)
    if annotations is not None:
        return annotations, "lrg", "lrg"

    annotations, reference_type = ncbi.fetch_annotations(reference_id, reference_type)
    if annotations is not None:
        return annotations, reference_type, "ncbi"

    annotations, reference_type = ensembl.fetch_annotations(
        reference_id, reference_type
    )
    if annotations is not None:
        return annotations, reference_type, "ensembl"
    return None, None, None


def fetch_sequence(reference_id, reference_source=None):
    """
    Sequence retrieval.

    :arg str reference_id: The reference ID.
    :arg str reference_source: The source from where to retrieve the sequence.

    :returns str: The sequence.
    """
    if reference_source is None:
        sequence = ncbi.fetch_sequence(reference_id)
        if sequence is None:
            lrg_annotations = lrg.fetch_lrg(reference_id)
            if lrg_annotations is not None:
                lrg_model = parser.parse(lrg_annotations, "lrg")
                if lrg_model is not None and lrg_model.get("sequence"):
                    sequence = lrg_model["sequence"]
        if sequence is None:
            sequence = ensembl.fetch_sequence(reference_id)
        return sequence
    else:
        if reference_source == "ncbi":
            return ncbi.fetch_sequence(reference_id)
        elif reference_source == "ensembl":
            return ensembl.fetch_sequence(reference_id)


def retrieve(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    parse=False,
    configuration_path=None,
):
    """
    Main retriever entry point. Identifies and calls the appropriate specific
    retriever methods.

    :arg str reference_id: The id of the reference.
    :arg str reference_source: The source of the reference, e.g., ncbi, ensembl.
    :arg str reference_type: The type of the reference: gff3, genbank, or lrg.
    :arg bool size_off: Flag for the maximum sequence length.
    :arg bool parse: Flag for parsing or not the reference.
    :arg bool parse: Flag for parsing or not the reference.
    :arg str configuration_path: Path towards configuration file.

    :return: The reference content and its type.
    """
    configuration.settings = configuration.setup_settings(configuration_path)

    annotations = None
    if reference_source is None and reference_type is None:
        annotations, reference_type, reference_source = fetch_annotations(
            reference_id, reference_type
        )
    elif reference_source is None and reference_type == "sequence":
        return fetch_sequence(reference_id)
    elif reference_source == "ncbi":
        if reference_type is None or reference_type == "gff3":
            annotations = ncbi.fetch_gff3(reference_id)
            reference_type = "gff3"
        elif reference_type == "genbank":
            annotations = ncbi.fetch_genbank(reference_id, not size_off)
            reference_type = "genbank"
        elif reference_type == "sequence":
            return fetch_sequence(reference_id, reference_source)
    elif reference_source == "ensembl":
        if reference_type is None or reference_type == "gff3":
            annotations = ensembl.fetch_gff(reference_id)
            reference_type = "gff3"
        elif reference_type == "json":
            annotations = ensembl.fetch_json(reference_id)
        elif reference_type == "sequence":
            return fetch_sequence(reference_id, reference_source)
    elif reference_source == "lrg":
        annotations = lrg.fetch_lrg(reference_id)

    if parse:
        if annotations is None:
            return
        model = parser.parse(annotations, reference_type, reference_source)
        if reference_type is "gff3":
            sequence = fetch_sequence(reference_id, reference_source)
            return {"model": model, "sequence": sequence, "source": reference_source}
        else:
            model.update({"source": reference_source})
            return model
    else:
        return annotations
