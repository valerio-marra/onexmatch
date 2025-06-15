from setuptools import setup, find_packages

setup(
    name='onexmatch',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'astropy',
    ],
    author='Valerio Marra',
    author_email='valerio.marra@me.com',
    description='Crossmatching tool for astronomical catalogs',
    python_requires='>=3.7',
)