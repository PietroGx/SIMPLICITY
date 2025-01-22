from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='simplicity',
    version='0.2.0',
    description='Stochastic sIMulation of SARS-CoV-2 sPread and evoLutIon aCcountIng for within-hosT dYnamics (SIMPLICITY)',
    long_description=readme,
    author='Pietro Gerletti, Jean-Baptiste Escudié',
    author_email='pietro.gerletti@pm.me',
    url='https://gitlab.com/combio/simplicity-public',
    license=license,
    packages=['simplicity'] 
    #packages=find_packages(exclude=('tests', 'docs'))
)
