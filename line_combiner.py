#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 14 08:01:36 2022

@author: kevin
"""

import numpy as np

with open('Galah_Blue_no_hfs.lin') as f:
    lines_all_no_hfs = f.readlines()
header_no_hfs=lines_all_no_hfs[:3]

lines_no_hfs=lines_all_no_hfs[3:]
line_names_no_hfs=lines_no_hfs[::4]

for value,y in enumerate(line_names_no_hfs):
    if len(y.split(',')[0])>20:
        ending_line=value
        break
line_names_no_hfs=line_names_no_hfs[:ending_line]
footer_no_hfs=lines_no_hfs[ending_line*4:]
lines_no_hfs=lines_no_hfs[:ending_line*4]

with open('Galah_Blue_hfs.lin') as f:
    lines_all_hfs = f.readlines()
header_hfs=lines_all_hfs[:3]

lines_hfs=lines_all_hfs[3:]
line_names_hfs=lines_hfs[::4]

for value,y in enumerate(line_names_hfs):
    if len(y.split(',')[0])>20:
        ending_line=value
        break
line_names_hfs=line_names_hfs[:ending_line]
footer_hfs=lines_hfs[ending_line*4:]
lines_hfs=lines_hfs[:ending_line*4]

new_lines=[x for x in line_names_hfs if not x in lines_no_hfs]