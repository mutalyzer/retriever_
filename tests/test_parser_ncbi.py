import pytest
import json
from pathlib import Path

from retriever.parsers.genbank import parse


def get_references():

    available_references = [
        'NM_078467.2',
        'NM_152263.2',
        'NM_152263.3',
        'NP_689476.2',
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
    assert json.loads(rmodel.loci_to_json_model()) == json.loads(model)


