#!/usr/bin/env python3

# ------------------------------------------------------------------
# Author: Libe Renteria Aizpurua
# Date: 2026-01-07 
#
# This script normalizes the lengths of the indices in a samplesheet CSV file by trimming the longer indices to match the shortest length.
# This is used because the downstream analysis requires all indices to be of the same length.
#
# ------------------------------------------------------------------

import csv
import sys
import argparse
from pathlib import Path

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--samplesheet', required=True, help='Path to samplesheet CSV')
args = parser.parse_args()

samplesheet = args.samplesheet

#load samplesheet
with open(samplesheet, 'r') as f:
    # Skip lines until we find [Data]
    for line in f:
        if line.strip() == '[Data]':
            break

    reader = csv.DictReader(f)
    rows = list(reader)
    fieldnames = reader.fieldnames

#get all index lengths
index_lengths = [len(row['index']) for row in rows]
max_length = max(index_lengths)
min_length = min(index_lengths)

modified_samples = []

#if we find a sample with a longer index, we will remove a nt
if max_length > min_length:
    for row in rows:
        if len(row['index']) == max_length:
            row['index'] = row['index'][:-1]
            modified_samples.append(row['Sample_ID'])
            print(f"Trimmed index for sample: {row['Sample_ID']}")

#now we write the modified samplesheet
samplesheet_path = Path(samplesheet)
output_name = f"{samplesheet_path.stem}_normalized.csv"
with open(output_name, 'w') as f:
    f.write('[Data]\n')
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

if modified_samples:
    with open('modified_samples.txt', 'w') as f:
        f.write('\n'.join(modified_samples))
        print(f"Modified {len(modified_samples)} samples")
else:
    print("No samples were modified")