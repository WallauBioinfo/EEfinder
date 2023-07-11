from setuptools import setup, find_packages

setup(
    name="eefinder",
    description="""
    Tool to analyze endogenous elements (EEs) based on similarity search
    and junctions of genomic regions.
    """,
    version="0.3.1",
    authors="Filipe Dezordi and Yago Dias",
    authors_emails='zimmer.filipe@gmail.com" and yag.dias@gmail',
    classifiers=[
        "Development Status :: Prototype",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click==8.1.3",
    ],
    entry_points={
        "console_scripts": [
            "eefinder = eefinder.scripts.main:main",
        ],
    },
)
