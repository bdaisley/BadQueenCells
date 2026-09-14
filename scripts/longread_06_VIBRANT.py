#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


def clean_value(key, value):
    """Remove surrounding quotes and normalize multiline qualifiers."""
    value = value.strip()

    if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
        value = value[1:-1]

    if key == "translation":
        value = re.sub(r"\s+", "", value)
    else:
        value = re.sub(r"\s+", " ", value).strip()

    return value


def looks_like_location(text):
    """Return True for GenBank-style feature locations."""
    compact = re.sub(r"\s+", "", text)

    if not compact:
        return False

    valid_start = (
        compact[0].isdigit()
        or compact.startswith("<")
        or compact.startswith(">")
        or compact.startswith("complement(")
        or compact.startswith("join(")
        or compact.startswith("order(")
    )

    return valid_start and bool(
        re.fullmatch(r"[A-Za-z0-9_<>\.\,\^\(\):]+", compact)
    )


if len(sys.argv) != 3:
    sys.exit(f"Usage: {Path(sys.argv[0]).name} input.gbk output.tsv")

input_gbk = Path(sys.argv[1])
output_tsv = Path(sys.argv[2])

rows = []
qualifier_names = set()

record_id = None
current_feature = None
current_qualifier = None
in_features = False


def save_feature():
    global current_feature, current_qualifier

    if current_feature is None:
        return

    if current_feature["feature_type"] != "CDS":
        current_feature = None
        current_qualifier = None
        return

    location = current_feature["location"]
    coordinates = [int(x) for x in re.findall(r"\d+", location)]

    if not coordinates:
        current_feature = None
        current_qualifier = None
        return

    row = {
        "record": current_feature["record"],
        "feature_type": current_feature["feature_type"],
        "location": location,
        "start": min(coordinates),
        "end": max(coordinates),
        "strand": -1 if location.startswith("complement(") else 1,
    }

    for key, values in current_feature["qualifiers"].items():
        cleaned = [clean_value(key, value) for value in values]
        row[key] = "; ".join(value for value in cleaned if value)
        qualifier_names.add(key)

    rows.append(row)

    current_feature = None
    current_qualifier = None


with input_gbk.open() as handle:
    for raw_line in handle:
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()

        # Record boundary
        match = re.match(r"^\s*LOCUS\s+(\S+)", line)
        if match:
            save_feature()
            record_id = match.group(1)
            in_features = False
            continue

        if re.match(r"^\s*FEATURES\b", line):
            save_feature()
            in_features = True
            continue

        if re.match(r"^\s*ORIGIN\b", line):
            save_feature()
            in_features = False
            continue

        if stripped == "//":
            save_feature()
            record_id = None
            in_features = False
            continue

        if not in_features or not stripped:
            continue

        # Feature lines, including tab-indented VIBRANT lines
        feature_match = re.match(
            r"^([A-Za-z][A-Za-z0-9_-]*)\s+(.+)$",
            stripped,
        )

        if feature_match:
            feature_type = feature_match.group(1)
            location = feature_match.group(2).strip()

            if looks_like_location(location):
                save_feature()

                current_feature = {
                    "record": record_id,
                    "feature_type": feature_type,
                    "location": re.sub(r"\s+", "", location),
                    "qualifiers": {},
                }

                current_qualifier = None
                continue

        if current_feature is None:
            continue

        # Qualifier line
        qualifier_match = re.match(
            r"^/([^=\s]+)(?:=(.*))?$",
            stripped,
        )

        if qualifier_match:
            key = qualifier_match.group(1)
            value = qualifier_match.group(2)

            if value is None:
                value = "true"

            current_feature["qualifiers"].setdefault(key, []).append(value)
            current_qualifier = key
            continue

        # Multiline qualifier continuation
        if current_qualifier is not None:
            current_feature["qualifiers"][current_qualifier][-1] += (
                " " + stripped
            )

save_feature()

base_columns = [
    "record",
    "feature_type",
    "location",
    "start",
    "end",
    "strand",
]

columns = base_columns + sorted(qualifier_names)

with output_tsv.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)

print(
    f"Wrote {len(rows)} CDS records from "
    f"{len(set(row['record'] for row in rows))} sequence records "
    f"with {len(columns)} columns"
)
