import io
from http.client import HTTPException
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

from .. import configuration
from ..util import make_request

Entrez.email = configuration.settings["EMAIL"]
Entrez.api_key = configuration.settings.get("NCBI_API_KEY")


class NoLinkError(Exception):
    """
    Raised when no transcript-protein (or vice-versa) link can be found.
    """

    pass


class NcbiConnectionError(Exception):
    """
    Raised when there is some NCBI connection error.
    """

    pass


class NoNcbiReference(Exception):
    """
    Raised when reference not found on NCBI.
    """

    pass


class ReferenceToLong(Exception):
    """
    Raised when the reference length exceeds maximum size.
    """

    pass


def fetch_ncbi_databases(reference_id):
    """
    Queries NCBI to identify in what databases a specific reference appears.

    Note:
    Whenever a reference without the most recent version is employed it seems
    that not all the databases are returned (including the most important ones
    for us, i.e., nuccore and protein). Hence, we always strip the reference
    and employ the accession only. Example:
    https://eutils.ncbi.nlm.nih.gov/gquery?term=NC_000001.10&retmode=xml
    versus:
    https://eutils.ncbi.nlm.nih.gov/gquery?term=NC_000001&retmode=xml
    Strange enough, but this works:
    https://www.ncbi.nlm.nih.gov/search/all/?term=NC_000001.10
    I was not able to find any additional parameter that would expand the
    search.

    :arg str reference_id: The id of the reference.
    :returns set: Set with NCBI databases.
    """
    if "." in reference_id:
        reference_id = reference_id.rsplit(".")[0]
    try:
        handle = Entrez.egquery(term=reference_id)
    except (IOError, HTTPError, HTTPException) as e:
        raise NcbiConnectionError

    result = Entrez.read(handle)
    databases = set()
    for item in result["eGQueryResult"]:
        if item["Status"].upper() == "OK" and int(item["Count"]) >= 1:
            databases.add(item["DbName"])
    return databases


def _fetch_reference_summary(reference_id):
    """
    Retrieves the reference summary if available on the NCBI.

    :arg str reference_id: The id of the reference.

    :returns dict:
    """
    try:
        databases = fetch_ncbi_databases(reference_id)
    except NcbiConnectionError as e:
        raise e

    if "nuccore" in databases:
        db = "nuccore"
    elif "protein" in databases:
        db = "protein"
    elif "nucest" in databases:
        db = "nucest"
        # Todo: Pay attention to the following:
        # https://ncbiinsights.ncbi.nlm.nih.gov/2018/07/30/upcoming-changes-est-gss-databases/
    else:
        raise NoNcbiReference

    try:
        handle = Entrez.esummary(db=db, id=reference_id)
    except (IOError, HTTPError, HTTPException) as e:
        print(
            "4, ERETR, Could not retrieve record length for {} from {}"
            "Entrez database. Error message {}:".format(reference_id, db, str(e))
        )
        raise NcbiConnectionError
    else:
        try:
            record = Entrez.read(handle)
        except RuntimeError:
            raise NoNcbiReference
        else:
            handle.close()

    return {
        "reference_id": record[0]["AccessionVersion"],
        "db": db,
        "length": int(record[0]["Length"]),
    }


def _fetch_link_from_ncbi(
    source_db,
    target_db,
    match_link_name,
    source_accession,
    source_version=None,
    match_version=True,
):
    """
    Retrieve a linked accession number from the NCBI.

    :arg str source_db: NCBI source database.
    :arg str target_db: NCBI target database.
    :arg function match_link_name: For each link found, this function is
      called with the link name (`str`) and it should return `True` iff the
      link is to be used.
    :arg str source_accession: Accession number for which we want to find a
      link (without version number).
    :arg int source_version: Optional version number for `source_accession`.
    :arg bool match_version: If `False`, the link does not have to match
      `source_version`.

    :raises NoLinkError: If no link could be retrieved from the NCBI.

    :returns: Tuple of `(target_accession, target_version)` representing the
      link target. If `source_version` is not specified or `match_version` is
      `False`, `target_version` can be `None`.
    :rtype: tuple(str, int)
    """

    # If we are currently strictly matching on version, we can try again if
    # no result is found. Otherwise, we just report failure.
    def fail_or_retry():
        if source_version is None or match_version:
            raise NoLinkError()
        return _fetch_link_from_ncbi(
            source_db,
            target_db,
            match_link_name,
            source_accession,
            source_version=None,
            match_version=False,
        )

    if source_version is None:
        source = source_accession
    else:
        source = "%s.%d" % (source_accession, source_version)

    # Find source record.
    try:
        handle = Entrez.esearch(db=source_db, term=source)
    except (IOError, HTTPException):
        # TODO: Log error.
        return fail_or_retry()

    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # TODO: Log error.
        return fail_or_retry()
    finally:
        handle.close()
    try:
        source_gi = result["IdList"][0]
    except IndexError:
        return fail_or_retry()

    # Find link from source record to target record.
    try:
        handle = Entrez.elink(dbfrom=source_db, db=target_db, id=source_gi)
    except (IOError, HTTPException):
        # TODO: Log error.
        return fail_or_retry()

    try:
        result = Entrez.read(handle)
    except Entrez.Parser.ValidationError:
        # TODO: Log error.
        return fail_or_retry()
    finally:
        handle.close()

    if not result[0]["LinkSetDb"]:
        return fail_or_retry()

    for link in result[0]["LinkSetDb"]:
        if match_link_name(link["LinkName"]):
            target_gi = link["Link"][0]["Id"]
            break
    else:
        return fail_or_retry()

    # Get target record.
    try:
        handle = Entrez.efetch(
            db=target_db, id=target_gi, rettype="acc", retmode="text"
        )
    except (IOError, HTTPException):
        # TODO: Log error.
        return fail_or_retry()

    target = handle.read().strip().split(".")
    handle.close()

    target_accession = target[0]
    target_version = int(target[1]) if source_version is not None else None
    return target_accession, target_version


def _fetch_link(
    source_db,
    target_db,
    match_link_name,
    source_accession,
    source_version=None,
    match_version=True,
):
    """
    The actual linking caller.
    """
    try:
        target_accession, target_version = _fetch_link_from_ncbi(
            source_db,
            target_db,
            match_link_name,
            source_accession,
            source_version=source_version,
            match_version=match_version,
        )
    except NoLinkError:
        raise

    return target_accession, target_version


def _transcript_to_protein(
    transcript_accession, transcript_version=None, match_version=True
):
    """
    Try to find the protein link to a transcript by using the NCBI Entrez API.

    :arg str transcript_accession: Accession number of the transcript for
      which we want to find the protein (without version number).
    :arg int transcript_version: Transcript version number. Please provide
      this if available, also if it does not need to match. This will enrich
      the cache.
    :arg bool match_version: If `False`, the link does not have to match
      `transcript_version`.

    :raises NoLinkError: If no link could be found.

    :returns: Tuple of `(protein_accession, protein_version)` representing the
      linked protein. If `transcript_version` is not specified or
      `match_version` is `False`, `protein_version` can be `None`.
    :rtype: tuple(str, int)
    """
    return _fetch_link(
        "nucleotide",
        "protein",
        lambda link: link in ("nuccore_protein", "nuccore_protein_cds"),
        transcript_accession,
        source_version=transcript_version,
        match_version=match_version,
    )


def _protein_to_transcript(protein_accession, protein_version=None, match_version=True):
    """
    Try to find the transcript link to a protein by using the NCBI Entrez API.

    :arg str protein_accession: Accession number of the protein for which we
      want to find the transcript (without version number).
    :arg int protein_version: Protein version number. Please provide this if
      available, also if it does not need to match. This will enrich the
      cache.
    :arg bool match_version: If `False`, the link does not have to match
      `protein_version`.

    :raises NoLinkError: If no link could be found.

    :returns: Tuple of `(transcript_accession, transcript_version)`
      representing the linked transcript. If `protein_version` is not
      specified or `match_version` is `False`, `transcript_version` can be
      `None`.
    :rtype: tuple(str, int)
    """
    return _fetch_link(
        "protein",
        "nucleotide",
        lambda link: link == "protein_nuccore_mrna",
        protein_accession,
        source_version=protein_version,
        match_version=match_version,
    )


def _link_transcript_to_protein_by_file(reference_id):
    """
    Try to find the link between a transcript and a protein (vice versa also)
    by downloading and using the corresponding genbank file.

    We only consider references starting with 'NM' (transcripts) and 'NP'
    (proteins), otherwise this is not reliable.

    For both NMs and NPs the feature with the answer is 'CDS'. Next, the
    qualifier for NMs is 'protein_id' (automatically in 'accession.version'
    format), while for NPs is 'coded_by' (it includes extra coordinates after
    the ':', .e.g, 'NM_012459.1:13..264').

    :arg str reference_id: The reference for which we try to get the link.
    :return: `accession[.version]` link reference.
    """
    if not (not reference_id.startswith("NP")) or (not reference_id.startswith("NM")):
        raise ValueError()
    content = fetch_genbank(reference_id)
    record = SeqIO.read(io.StringIO(content), "genbank")
    for feature in record.features:
        if feature.type == "CDS":
            if feature.qualifiers.get("protein_id"):
                protein_id = feature.qualifiers["protein_id"][0]
            elif feature.qualifiers.get("coded_by"):
                protein_id = feature.qualifiers["coded_by"][0]
                if ":" in protein_id:
                    protein_id = protein_id.rsplit(":")[0]
            else:
                raise NoLinkError()
            return _decompose_reference(protein_id)

    raise NoLinkError()


def link_reference(reference_id):
    """
    Make the appropriate calls to return the protein/transcript link.

    :arg str reference_id: The reference for which we try to get the link.
    :return: `accession[.version]` link reference.
    :rtype: str
    """
    accession, version = _decompose_reference(reference_id)
    if version is None:
        match_version = False
    else:
        match_version = True

    try:
        link_accession, link_version = _protein_to_transcript(
            accession, version, match_version
        )
    except NoLinkError:
        pass
    else:
        if link_accession:
            return _compose_reference(link_accession, link_version), "api"

    try:
        link_accession, link_version = _transcript_to_protein(
            accession, version, match_version
        )
    except NoLinkError:
        pass
    else:
        return _compose_reference(link_accession, link_version), "api"

    if version:
        try:
            link_accession, link_version = _link_transcript_to_protein_by_file(
                reference_id
            )
        except NoLinkError:
            return None, None
        except ValueError:
            return None, None
        else:
            return _compose_reference(link_accession, link_version), "file"

    return None, None


def _decompose_reference(reference_id):
    """
    Get the accession and the version of a reference. The version is None
    if it is not present.
    """
    if "." in reference_id:
        accession, version = reference_id.rsplit(".", 1)
        version = int(version)
    else:
        accession = reference_id
        version = None
    return accession, version


def _compose_reference(accession, version=None):
    """
    Create the accession[.version] of a reference.
    """
    if accession is None:
        return None
    if version is None:
        return accession
    else:
        return "{}.{}".format(accession, version)


def fetch_genbank(reference_id, size_on=True):
    """
    Retrieve a genbank reference from the NCBI.

    :arg str reference_id: The id of the reference.
    :arg bool size_on: Consider or not the maximum sequence length.
    :returns: Reference content.
    :rtype: str
    """
    try:
        reference_summary = _fetch_reference_summary(reference_id)
    except NoNcbiReference:
        return None
    except NcbiConnectionError:
        raise NcbiConnectionError

    if (
        size_on
        and reference_summary["length"] > configuration.settings["MAX_FILE_SIZE"]
    ):
        raise ReferenceToLong
    try:
        handle = Entrez.efetch(
            db=reference_summary["db"],
            id=reference_id,
            rettype="gbwithparts",
            retmode="text",
        )
    except (IOError, HTTPError, HTTPException):
        raise NcbiConnectionError
    else:
        raw_data = handle.read()
        handle.close()
        return raw_data


def fetch_sequence(reference_id):
    """
    Retrieve the sequence of corresponding reference ID.

    :arg str reference_id: The reference ID.

    :returns str: The sequence.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="fasta")
    except HTTPError:
        return
    records = []
    for record in SeqIO.parse(handle, "fasta"):
        records.append({"seq": str(record.seq), "description": record.description})
    handle.close()
    if len(records) == 1:
        return records[0]


def fetch_gff3(reference_id):
    """
    Retrieve the gff3 for the corresponding reference ID.

    :arg str reference_id: The reference ID.

    :returns str: gff3 content.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
    params = {"db": "nuccore", "report": "gff3", "id": reference_id}
    return make_request(url, params)


def fetch_annotations(reference_id, reference_type="gff3", size_on=True):
    """
    Fetch the raw annotation for the corresponding reference ID.

    :arg str reference_id: The reference ID.
    :arg str reference_type: The reference type ("gff3" - default, or "genbank").
    :arg bool size_on: Consider maximum file size.

    :returns tuple: raw annotations, type ("gff3" or "genbank")
    """
    if reference_type in [None, "gff3"]:
        return fetch_gff3(reference_id), "gff3"
    if reference_type == "genbank":
        return fetch_genbank(reference_id, size_on), "genbank"
    return None, None
