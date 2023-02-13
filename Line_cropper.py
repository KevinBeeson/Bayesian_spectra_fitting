#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:46:20 2021

@author: kevin
"""

import csv
import numpy as np


list_stars = list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))
original_linelist=list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))