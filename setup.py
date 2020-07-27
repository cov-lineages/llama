from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from llama import __version__, _program

setup(name='llama',
      version=__version__,
      packages=find_packages(),
      scripts=["llama/scripts/Snakefile",
      "llama/scripts/assess_putative_lineage.smk",
      "llama/scripts/check_metadata.py",
      "llama/scripts/find_closest_in_db.smk",
      "llama/scripts/parse_paf.py",
      "llama/scripts/process_local_trees.smk",
      "llama/scripts/just_collapse_trees.smk",
      "llama/scripts/make_report.py",
      "llama/scripts/make_tree_figures.py",
      "llama/scripts/data_parsing.py",
      "llama/scripts/baltic.py",
      "llama/scripts/report_template.pmd"],
      package_data={"llama":["data/reference.fasta",
                              "data/footer.png"]},
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            "pweave>=0.30.3",
            "matplotlib>=3.2.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "scipy>=1.4.1",
            "numpy>=1.13.3"
        ],
      description='Local Lineage and Monophyly Assessment',
      url='github.com/cov-lineages/llama',
      author='Aine OToole & Rambaut Lab',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = llama.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
