# hkbin
# Copyright (C) 2020-2021  Sebastiaan L. Zoutendijk
# Based on pyGravSphere
# Copyright (C) 2019-2020  Anna Genina, Justin Read
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys

galaxies = sys.argv[1:]

codedir = str(os.environ["GravSpherePath"])
if codedir[-1] != '/':
    codedir  = codedir + '/'
sys.path.append(codedir)

import hkgal_input as gal_input

for i, galaxy in enumerate(galaxies):
    print("Preprocessing galaxy {:d}/{:d}: {:s}".format(i+1, len(galaxies), galaxy))
    gal_input.galaxy_data(galaxy, "GalaxyData/", "KinPhotDat/")
