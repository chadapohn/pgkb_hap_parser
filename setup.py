from setuptools import setup, find_packages

setup(
    name = 'pgkb_hap_parser',
    version = '0.0.1',
    packages = find_packages(),
    entry_points = {
        'console_scripts': [
            'pgkb_hap_parser = main:parser',
        ],
    },
)