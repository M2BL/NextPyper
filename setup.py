import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "nextpyper",
            "nextpyper.VERSION",
        )
    ) as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Other/Proprietary License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="nextpyper",
    packages=find_packages(),
    url="https://git.sorbus.ibot.cas.cz/m2b_ibot/nextpyper",
    python_requires=">=3.12",
    description="Recovery of homologous genes from targeted sequence capture data for higher ploidy samples",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Simón Villanueva Corrales",
    author_email="simon.corrales@ibot.cas.cz",
    data_files=get_data_files(),
    py_modules=["nextpyper"],
    install_requires=[
        "snaketool-utils>=0.0.4",
        "snakemake>=9.3.1",
        "pyyaml>=6.0",
        "Click>=8.1.3",
        "pandas>=2.2",
        "biopython>=1.83",
        "numpy>=1.26",
        "scikit-learn>=1.5",
        "scipy>=1.14",
        "polars>=1.9",
        "intervaltree>=3.1",
        "kmedoids>=0.5.3.1",
        "more-itertools>=10.7",
    ],
    entry_points={"console_scripts": ["nextpyper=nextpyper.__main__:main"]},
    include_package_data=True,
)
