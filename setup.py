"""
Setup script for GCMCbyGULP
===========================

"""

from setuptools import setup, find_packages

setup(
    name='GCMCbyGULP',
    version='0.1.0',
    description='Utilities for performing GCMC simulations by GULP',
    author='Tschijnmo TSCHAU',
    author_email='tschijnmotschau@gmail.com',
    url='https://github.com/tschijnmo/GCMCbyGULP',
    license='MIT',

    packages=find_packages(),
    install_requires=[
        'CoolProp',
        'periodic',
        'pystache',
        'numpy',
        'scipy',
        'PyYAML',
        ],
    entry_points={
        'console_scripts': [
            'GCMCbyGULP-gen = GCMCbyGULP.geninp:gen_main',
            'GCMCbyGULP-read = GCMCbyGULP.readres:read_main',
            ],
        },
    )
