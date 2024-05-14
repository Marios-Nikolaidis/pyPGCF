from setuptools import setup

__version__ = "0.1.4"

setup(
    name="pypgcf",
    version=__version__,
    author="Marios Nikolaidis",
    author_email="marionik23@gmail.com",
    url="https://github.com/m-nikolaidis/pyPGCF",
    license="GPL",
    description="pyPGCF: PhyloGenomic, Core and Fingerprint analyses package",
    packages=["pypgcf"],
    entry_points={
        "console_scripts": [
            "pyPGCF = pypgcf.cli:main",
        ],
    },
)
