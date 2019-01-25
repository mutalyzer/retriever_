import pytest
import json
from pathlib import Path

from retriever.parsers.genbank import parse


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def get_references():

    available_references = [
        'NM_078467.2',
        'NM_152263.2',
        'NM_152263.3',
        'NP_689476.2',
        'NG_012337.1',
        'NR_002196.2',
        'L41870.1'
    ]

    references = []

    for reference in available_references:

        path_gb = Path(Path(__file__).parent) / 'data' / str(reference + '.gb')
        with path_gb.open() as f:
            gb = f.read()

        path_json = Path(Path(__file__).parent) / 'data' / str(reference + '.json')
        with path_json.open() as f:
            json = f.read()

        references.append((reference, gb, json))

    return references


@pytest.mark.parametrize(
    'reference, content, model', get_references()
)
def test_model(reference, content, model):
    rmodel = parse(content)
    assert ordered(rmodel.to_dict()) == ordered(json.loads(model))
