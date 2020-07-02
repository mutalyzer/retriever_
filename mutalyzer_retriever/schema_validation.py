from schema import And, Optional, Schema

point = Schema(
    {
        "type": "point",
        "position": And(int, lambda n: n >= 0),
        Optional("uncertain"): And(bool, True),
    }
)

location = Schema({"type": str, "start": point, "end": point, Optional("strand"): int})

feature = Schema(
    {
        "type": str,
        "id": str,
        Optional("location"): location,
        Optional("features"): list,
        Optional("qualifiers"): dict,
    }
)


def validate(current):
    feature.validate(current)
    if current.get("features"):
        for sub_feature in current["features"]:
            validate(sub_feature)
