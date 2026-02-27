# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='simplicity',
    version='2.2.0',
    description='Stochastic sIMulation of SARS-CoV-2 sPread and evoLutIon aCcountIng for within-hosT dYnamics (SIMPLICITY)',
    long_description=readme,
    author='Pietro Gerletti, Jean-Baptiste Escudié',
    author_email='simplicity.founder544@passmail.com',
    url='https://github.com/PietroGx/SIMPLICITY',
    license=license,
    packages=['simplicity'] 
    #packages=find_packages(exclude=('tests', 'docs'))
)
