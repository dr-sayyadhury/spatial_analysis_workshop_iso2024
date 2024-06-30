from setuptools import setup, find_packages

setup(
    name='spatial_pre_processing_functions',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'seaborn',
        'pandas',
        'pyarrow',
        'fastparquet',
    ],
)

