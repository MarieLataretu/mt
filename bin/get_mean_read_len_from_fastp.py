#!/usr/bin/env python3

import sys

import json

json_report = sys.argv[1]
mode = sys.argv[2]

with open(json_report) as json_report:
    report = json.load(json_report)

if mode == 'single':
    mean_read_length = report['summary']['after_filtering']['read1_mean_length']
elif mode == 'paired':
    mean_read_length = (float(report['summary']['after_filtering']['read1_mean_length']) + float(report['summary']['after_filtering']['read2_mean_length'])) / 2
else:
    sys.exit(f"Invalid read type: {mode}")

print(mean_read_length)