#!/usr/bin/env python3

#!/usr/bin/env python3

# ------------------------------------------------------------------
# Author: Libe Renteria Aizpurua
# Date: 2026-01-07 
#
# This script creates a demultiplexing samplesheet for the fqtk demultiplexing tool.
# It takes the forward sequence of the second index to create the full barcode, because the selected sequencing instrument requires it.
#
# ------------------------------------------------------------------

import csv
import sys
import argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--samplesheet', required=True, help='Path to samplesheet CSV')
args = parser.parse_args()

samplesheet = args.samplesheet

# Shared parser accepts Illumina [Data] sections or a plain CSV header.
from validate_inputs import read_samplesheet
rows = read_samplesheet(samplesheet)
fieldnames = list(rows[0])

#create demux_sheet
with open('DemuxSampleSheet.tsv', 'w') as f:
    f.write('sample_id\tbarcode\n')
    for row in rows:
        sample_id = row['Sample_ID']
        barcode = row['index'] + row['index2']
        f.write(f'{sample_id}\t{barcode}\n')