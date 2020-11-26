import json

from ..request import Http400, RequestErrors, request


def fetch_json(feature_id, timeout=1):
    url = "https://rest.ensembl.org/lookup/id/{}".format(feature_id)
    params = {"feature": ["gene", "transcript", "cds"], "expand": 1}
    headers = {"Content-Type": "application/json"}
    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError
        else:
            raise e
    else:
        return response


def fetch_sequence_details(feature_id, timeout=1):
    url = "https://rest.ensembl.org/lookup/id/{}".format(feature_id)
    headers = {"Content-Type": "application/json"}
    response = json.loads(request(url, headers=headers, timeout=timeout))
    return (
        response["start"],
        response["end"],
        response["species"],
        response["seq_region_name"],
    )


def fetch_sequence(feature_id, timeout=1):
    start, end, species, seq_region_name = fetch_sequence_details(feature_id, timeout)
    url = "https://rest.ensembl.org/sequence/region/{}/{}:{}..{}".format(
        species, seq_region_name, start, end
    )
    headers = {"Content-Type": "application/json"}
    return json.loads(request(url, headers=headers, timeout=timeout))


def fetch_fasta(feature_id, timeout=1):
    url = "https://rest.ensembl.org/sequence/id/{}".format(feature_id)
    params = {"format": "fasta", "type": "genomic"}
    headers = {"Content-Type": "text/x-fasta"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError
        else:
            raise e
    else:
        return response


def fetch_gff3(feature_id, timeout=1):
    url = "https://rest.ensembl.org/overlap/id/{}".format(feature_id)
    params = {"feature": ["gene", "transcript", "cds", "exon"]}
    headers = {"Content-Type": "text/x-gff3"}

    try:
        response = request(url, params, headers, timeout=timeout)
    except RequestErrors:
        raise ConnectionError
    except Http400 as e:
        response_json = e.response.json()
        if response_json and response_json.get("error") == "ID '{}' not found".format(
            feature_id
        ):
            raise NameError
        else:
            raise e
    else:
        return response


def fetch(reference_id, reference_type, timeout=1):
    if reference_type in [None, "gff3"]:
        return fetch_gff3(reference_id, timeout), "gff3"
    elif reference_type == "fasta":
        return fetch_fasta(reference_id, timeout), "fasta"
    elif reference_type == "json":
        return fetch_json(reference_id, timeout), "json"
    elif reference_type == "genbank":
        return None, "genbank"

    raise ValueError(
        "Ensembl fetch does not support '{}' reference type.".format(reference_type)
    )
