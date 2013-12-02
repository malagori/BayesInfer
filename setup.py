#!/usr/bin/env python

from setuptools import setup
setup(
    name='BayesInfer',
    version='0.1dev',
    author='Mehmood Alam Khan',
    author_email='malagori@kth.se',
    url='https://github.com/malagori/BayesInfer',
    packages=['bayesInfer',],
    scripts = ['fuv','scripts/dfToBeneFormat.py'],
    license='GPLv3',
    long_description=open('README.md').read(),
    install_requires= ['numpy >=1.6.1'],
)