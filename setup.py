from setuptools import setup, find_packages
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from llama import __version__, _program

setup(name='llama',
      version=__version__,
      packages=find_packages(),
      scripts=["llama/scripts/Snakefile",
      "llama/scripts/parse_paf.py",
      "llama/scripts/find_closest_in_tree.smk",
      "llama/scripts/assess_input_file.smk",
      "llama/scripts/process_collapsed_trees.smk",
      "llama/scripts/process_catchment_trees.smk",
      "llama/scripts/just_collapse_trees.smk"],
      package_data={"llama":["data/reference.fasta"},
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "matplotlib>=3.2.1",
            "pweave>=0.30.3",
            "scipy>=1.4.1",
            "numpy>=1.13.3"
        ],
      description='Local Lineage and Monophyly Assessment',
      url='https:cov-lineages.org/llama',
      author='Aine OToole & Andrew Rambaut',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = llama.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
