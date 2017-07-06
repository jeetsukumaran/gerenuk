#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="gerenuk",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["gerenuk", "test"],
    # scripts=["bin/gerenuk.py",],
    url="http://pypi.python.org/pypi/gerenuk/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    # install_requires=[ ],
)