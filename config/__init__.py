#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Configurations module
'''
# ~~~~~ SETUP ~~~~~~ #
import yaml
import os

scriptdir = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(scriptdir, "sns.yml"), "r") as f:
    sns = yaml.load(f)
