import json
from urllib.error import HTTPError

from ..util import make_request


def fetch_json(feature_id):
    url = "https://rest.ensembl.org/lookup/id/{}".format(feature_id)
    params = {"feature": ["gene", "transcript", "cds"], "expand": 1}
    headers = {"Content-Type": "application/json"}
    return make_request(url, params, headers)


def fetch_gff(feature_id):
    url = "https://rest.ensembl.org/overlap/id/{}".format(feature_id)
    params = {"feature": ["gene", "transcript", "cds", "exon"]}
    headers = {"Content-Type": "text/x-gff3"}
    try:
        return make_request(url, params, headers)
    except HTTPError as err:
        print("HTTP error")


def fetch_sequence_details(feature_id):
    url = "https://rest.ensembl.org/lookup/id/{}".format(feature_id)
    headers = {"Content-Type": "application/json"}
    response = json.loads(make_request(url, headers=headers))
    return (
        response["start"],
        response["end"],
        response["species"],
        response["seq_region_name"],
    )


def fetch_sequence(feature_id):
    start, end, species, seq_region_name = fetch_sequence_details(feature_id)
    url = "https://rest.ensembl.org/sequence/region/{}/{}:{}..{}".format(
        species, seq_region_name, start, end
    )
    headers = {"Content-Type": "application/json"}
    return json.loads(make_request(url, headers=headers))


def fetch_annotations(reference_id, reference_type):
    if reference_type in [None, "gff3"]:
        return fetch_gff(reference_id), "gff3"
    if reference_type == "json":
        return fetch_json(reference_id), "json"
    if reference_type == "genbank":
        return None, "genbank"
    return None, None
