from http.client import HTTPException
from http.client import HTTPException
from urllib.error import HTTPError

from Bio import Entrez

# Maximum size for uploaded and downloaded files (in bytes).
MAX_FILE_SIZE = 10 * 1048576  # 10 MB

# This email address is used in contact information on the website and sent
# with NCBI Entrez calls.
EMAIL = 'mutalyzer@lumc.nl'

Entrez.email = EMAIL


class NoLinkError(Exception):
    """
    Raised when no transcript-protein (or vice-versa) link can be found.
    """
    pass


def length(reference_id, db='nuccore'):
    """
    Retrieves the sequence length.

    :arg str reference_id: The id of the reference.
    :returns: Sequence length.
    :rtype: int
    """
    try:
        net_handle = Entrez.esummary(db=db, id=reference_id)
        record = Entrez.read(net_handle)
        net_handle.close()
    except (IOError, HTTPError, HTTPException) as e:
        print('-1, INFO, Error connecting to Entrez {} database: {}.'
              .format(db, e))
        print('4, ERETR, Could not retrieve record length for {}.'
              .format(reference_id))
        return None
    except RuntimeError as e:
        # TODO: Extract the db name from the error:
        # Otherdb uid="131889391" db="protein" term="131889391"
        return length(reference_id, db='protein')
    return int(record[0]['Length'])


def fetch_ncbi(reference_id, size_on=True):
    """
    Retrieves a genbank reference from the NCBI.

    :arg str reference_id: The id of the reference.
    :arg bool size_on: Consider or not the maximum sequence length.
    :returns: Reference content.
    :rtype: str
    """
    if size_on:
        if length(reference_id) is None:
            print('Retrieval stopped since there is no record length.')
            return
        if length(reference_id) > MAX_FILE_SIZE:
            print('4, ERETR, Could not retrieve {} (exceeds maximum file '
                  'size of {} megabytes).'
                  .format(reference_id, MAX_FILE_SIZE // 1048576))
            return None
    try:
        net_handle = Entrez.efetch(
            db='nuccore', id=reference_id, rettype='gbwithparts',
            retmode='text')
        raw_data = net_handle.read()
        net_handle.close()
    except (IOError, HTTPError, HTTPException) as e:
        print('-1, INFO, Error connecting to Entrez nuccore database: {}.'
              .format(e))
        print('4, ERETR, Could not retrieve record {}.'
              .format(reference_id))
        return None

    return raw_data


def _get_link_from_ncbi(source_db, target_db, match_link_name,
                        source_accession, source_version=None,
                        match_version=True):
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
        return _get_link_from_ncbi(source_db, target_db, match_link_name,
                                   source_accession, source_version=None,
                                   match_version=False)
    if source_version is None:
        source = source_accession
    else:
        source = '%s.%d' % (source_accession, source_version)

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
        source_gi = result['IdList'][0]
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

    if not result[0]['LinkSetDb']:
        return fail_or_retry()

    for link in result[0]['LinkSetDb']:
        if match_link_name(link['LinkName']):
            target_gi = link['Link'][0]['Id']
            break
    else:
        return fail_or_retry()

    # Get target record.
    try:
        handle = Entrez.efetch(
            db=target_db, id=target_gi, rettype='acc', retmode='text')
    except (IOError, HTTPException):
        # TODO: Log error.
        return fail_or_retry()

    target = handle.read().strip().split('.')
    handle.close()

    target_accession = target[0]
    target_version = int(target[1]) if source_version is not None else None
    return target_accession, target_version


def _get_link(source_db, target_db, match_link_name,
              source_accession, source_version=None, match_version=True):
    """
    The actual linking caller.
    """
    try:
        target_accession, target_version = _get_link_from_ncbi(
            source_db, target_db, match_link_name, source_accession,
            source_version=source_version, match_version=match_version)
    except NoLinkError:
        raise

    return target_accession, target_version


def transcript_to_protein(transcript_accession, transcript_version=None,
                          match_version=True):
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
    return _get_link(
        'nucleotide', 'protein',
        lambda link: link in ('nuccore_protein', 'nuccore_protein_cds'),
        transcript_accession, source_version=transcript_version,
        match_version=match_version)


def protein_to_transcript(protein_accession, protein_version=None,
                          match_version=True):
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
    return _get_link(
        'protein', 'nucleotide', lambda link: link == 'protein_nuccore_mrna',
        protein_accession, source_version=protein_version,
        match_version=match_version)


def link_reference(reference_id):
    """
    Make the appropriate calls to return the protein/transcript link.

    :arg str reference_id: The reference for which we try to get the link.
    :return: `accession[.version]` link reference.
    :rtype: str
    """
    if '.' in reference_id:
        accession, version = reference_id.rsplit('.', 1)
        version = int(version)
        match_version = True
    else:
        accession = reference_id
        version = None
        match_version = False

    try:
        link_accession, link_version = protein_to_transcript(
            accession, version, match_version)
    except NoLinkError:
        pass
    else:
        if link_accession:
            return compose_reference(link_accession, link_version)

    try:
        link_accession, link_version = transcript_to_protein(
            accession, version, match_version)
    except NoLinkError:
        pass
    else:
        return compose_reference(link_accession, link_version)

    return None


def compose_reference(accession, version=None):
    """
    Get the accession[.version] of a reference.
    """
    print(accession, version)
    if accession is None:
        return None
    if version is None:
        return accession
    else:
        return '{}.{}'.format(accession, version)
