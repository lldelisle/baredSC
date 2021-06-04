# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup


def get_version():
  try:
    f = open(os.path.join("baredSC", "_version.py"))
  except EnvironmentError:
    return None
  for line in f.readlines():
    mo = re.match("__version__ = '([^']+)'", line)
    if mo:
      ver = mo.group(1)
      return ver
  return None


# From numpy:
def get_docs_url():
  # For releases, this URL ends up on pypi.
  # By pinning the version, users looking at old PyPI releases can get
  # to the associated docs easily.
  return f"https://baredsc.readthedocs.io/en/v{get_version()}"


install_requires_py = ["numpy >=1.16",
                       "matplotlib >=3.1.1",
                       "pandas >=0.25.0",
                       "scipy >=1.3.0",
                       "corner >=2.0.0",
                       "samsam >=0.1.2",
                       "anndata >=0.7"
                       ]


setup(
    name='baredSC',
    version=get_version(),
    author='Lucille Lopez-Delisle, Jean-Baptiste Delisle',
    author_email='lucille.delisle@epfl.ch',
    packages=['baredSC'],
    scripts=['bin/baredSC_1d', 'bin/baredSC_2d',
             'bin/combineMultipleModels_1d',
             'bin/combineMultipleModels_2d'],
    url='https://github.com/lldelisle/baredSC',
    project_urls={
        "Bug Tracker": "https://github.com/lldelisle/baredSC/issues",
        "Documentation": get_docs_url()
    },
    license='GPLv3',
    description='baredSC: Bayesian Approach to Retreive Expression Distribution of Single Cell',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    classifiers=[
                'Intended Audience :: Science/Research',
                'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=install_requires_py,
    zip_safe=False,
    python_requires='>=3.7.*, <4'
)
