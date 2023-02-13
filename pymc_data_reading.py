#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:37:09 2022

@author: kevin
"""

import arviz as az



data=az.from_netcdf('4_cores_1600')
az.summary(data)
az.plot_trace(data)