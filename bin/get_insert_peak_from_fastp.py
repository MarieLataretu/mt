#!/usr/bin/env python3

import sys

import json

json_report = sys.argv[1]

with open(json_report) as json_report:
    report = json.load(json_report)

insert_size = report['insert_size']['peak']

print(insert_size)