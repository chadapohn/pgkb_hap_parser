#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = "pgkb_hap_parser",
    version = "0.1.0",
    packages = find_packages(),
    entry_points = {
        "console_scripts": [
            "pgkb_hap_parser = src.pgkb_hap_parser.parser:parse",
        ]
    }
)
