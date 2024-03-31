#!/usr/bin/env python

from distutils.core import setup

setup(
    name='PySTFP',
    version='1.0.0',
    url='https://github.com/juliankappler/PySTFP',
    author='Julian Kappler',
    author_email='jkappler@posteo.de',
    description='Python module for the 1D Fokker-Planck short-time propagator',
    packages=['PySTFP',
    'PySTFP.sympy_definitions',
    ],
)
