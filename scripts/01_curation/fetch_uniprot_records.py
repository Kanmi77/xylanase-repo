#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import argparse
import csv
import re
import time

import yaml


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

UNIPROT_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "organism_name",
    "lineage",
    "length",
    "sequence",
    "xref_cazy",
    "xref_pdb",
    "reviewed",
]


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def fetch_page(url):
    request = Request(
        url,
        headers={"User-Agent": "enzyme-thermostability-workflow"},
    )

    with urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8")
        link_header = response.headers.get("Link", "")

    return text, link_header


def get_next_url(link_header):
    """
    Extract the UniProt pagination URL from the Link header.

    Do not split the header by comma, because the fields query parameter
    can also contain commas.
    """
    match = re.search(r"<([^>]+)>;\s*rel=\"next\"", link_header)

    if match:
        return match.group(1)

    return None


def parse_tsv(text):
    lines = [line for line in text.splitlines() if line.strip()]

    if not lines:
        return [], []

    reader = csv.DictReader(lines, delimiter="\t")
    return reader.fieldnames, list(reader)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch UniProt records for the configured enzyme query."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    query = config["enzyme_search"]["uniprot_query"]

    parameters = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(UNIPROT_FIELDS),
        "size": 500,
    }

    url = f"{UNIPROT_SEARCH_URL}?{urlencode(parameters)}"
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    fieldnames = None
    page_count = 0

    while url:
        page_count += 1
        print(f"Fetching page {page_count}: {url}")

        text, link_header = fetch_page(url)
        page_fields, rows = parse_tsv(text)

        if fieldnames is None:
            fieldnames = page_fields

        all_rows.extend(rows)
        url = get_next_url(link_header)

        if url:
            time.sleep(0.2)

    if not fieldnames:
        raise SystemExit("No UniProt records were returned.")

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Fetched pages: {page_count}")
    print(f"Fetched records: {len(all_rows)}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
