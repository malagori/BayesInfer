#!/usr/bin/env python

from setuptools import setup
setup(
    name='malagori',
    version='0.1dev',
    author='Mehmood Alam Khan',
    author_email='malagori@kth.se',
    url='https://github.com/malagori/BayesInfer',
    packages=['bayesInfer',],
    scripts = ['fuv','scripts/calculateBDeuScore.py'],
    license='GPLv3',
    long_description=open('README.md').read(),
)