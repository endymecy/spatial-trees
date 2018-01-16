#!/usr/bin/env Python
# -*- coding=utf-8 -*-

from distutils.core import setup
from setuptools import find_packages

setup(
    name = 'spatial-trees',
    version = '0.0.1',
    description = 'Implement of some spatial tree',
    long_description=open("README.md").read(),
    author = 'endymecy',
    author_email = 'endymecy@sina.cn',
    url = 'https://github.com/endymecy/spatial-trees',
    license="Apache License 2.0",
    keywords='spatial trees',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
    install_requires=['numpy','scipy'],
)