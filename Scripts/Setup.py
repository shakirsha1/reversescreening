#!/usr/bin/env python3
"""
PharmacoNet Setup
"""

from setuptools import setup, find_packages

setup(
    name="pmnet",
    version="2.2.0",
    description="PharmacoNet: Deep Learning-based Pharmacophore Modeling",
    author="PharmacoNet Team",
    author_email="your-email@example.com",
    url="https://github.com/your-username/pharmaconet",
    packages=find_packages(),
    install_requires=[
        "torch>=2.0.0",
        "rdkit>=2022.9.1",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "matplotlib>=3.6.0",
        "seaborn>=0.12.0",
        "scikit-learn>=1.1.0",
        "biopython>=1.80",
    ],
    python_requires=">=3.11",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
    ],
)
