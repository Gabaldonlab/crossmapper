#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 13:19:26 2018

@author: hrant
"""

import sys
import subprocess

cmd_norm = "grep '>' ../test1.fasta"
cmd_error = "grep '>' no_file"



subprocess.run(cmd_error, shell = True, capture_output=True)


subprocess.run(cmd_norm, shell = True)

