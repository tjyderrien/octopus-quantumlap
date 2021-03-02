#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for spp_extended_theory.
    Use setup.cfg to configure your project.
"""
import sys

from pkg_resources import require, VersionConflict
from setuptools import setup

try:
    require("setuptools>=38.3")
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)

try:
    require("spp_extended_theory")
except VersionConflict:
    print("Error: version of spp_extended_theory is ill-defined!")
    sys.exit(1)

if __name__ == "__main__":
    setup()
