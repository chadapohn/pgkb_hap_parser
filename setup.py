from setuptools import setup, find_packages

setup(
    name = 'pgkb_hap_parser',
    version = '0.0.7',
    packages = find_packages(),
    entry_points = {
        'console_scripts': [
            'pgkb_hap_parser = pgkb_hap_parser.parser:parse',
        ],
    },
)