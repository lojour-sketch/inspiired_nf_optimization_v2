#!/usr/bin/env python3

import csv
import sys
import argparse

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

#create demux_sheet
with open('DemuxSampleSheet.tsv', 'w') as f:
    f.write('sample_id\tbarcode\n')
    for row in rows:
        sample_id = row['Sample_ID']
        barcode = row['index'] + row['index2']
        f.write(f'{sample_id}\t{barcode}\n')