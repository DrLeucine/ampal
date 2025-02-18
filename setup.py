"""Setup script for the AMPAL framework."""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize


setup(
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    setup_requires=[
        "Cython",
    ],
    ext_modules=cythonize(
        [
            Extension("ampal.geometry", ["src/ampal/geometry.pyx"]),
        ]
    ),
    zip_safe=False,
)
