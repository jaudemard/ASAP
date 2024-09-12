"""
Setup file to run asap.

In a terminal and from the ASAP repository (where setup.py is located), run:

pip install .
"""
import os
import sys
from setuptools import setup, find_packages

base_path = os.path.dirname(__file__)
sys.path.insert(0, base_path)

setup(
    name="asap",
    version="24.09",
    author="Juliette AUDEMARD",
    author_email="juliette.aud@live.fr",
    description="Compute solvent accessible surface area (SASA) using Shrake-Rupley method.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=['asap'],
    package_dir={'asap': 'src/asap'},
    entry_points={
        'console_scripts': [
            'asap = asap.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.7',
)
