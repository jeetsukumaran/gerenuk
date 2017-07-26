#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="gerenuk",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["gerenuk"],
    scripts=[
        "bin/gerenuk-simulate.py",
        "bin/gerenuk-reject.py",
        ],
    url="http://pypi.python.org/pypi/gerenuk/",
    test_suite = "gerenuk.test",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    # install_requires=[ ],
)
