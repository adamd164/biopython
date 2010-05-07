#!/usr/bin/env python
#
# CoordBuilder.py
"""
Implementation of the NeRF algorithm
for constructing cartesian from torsion space coords.

Approach is described in Parsons, Strauss Torsion Space .. JCC 2005
"""
#TODO: complete reference

#TODO: Code ready for BioPython. KR 2010/05/07

#TODO: re-license

__author__ = "Kristian Rother, Magdalena Musielak, Tomasz Puton"
__copyright__ = "Copyright 2008, The Moderna Project"
__license__ = "GPL"
__credits__ = ["Janusz Bujnicki"]
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from numpy import array, ndarray
from Bio.PDB.Vector import Vector
import math

def get_ref_matrix(vec1, vec2, vec3):
    """
    Returns a matrix to transform reference coordinates into real coordinates.
    """
    vec23 = (vec3 - vec2) / (vec2 - vec3).norm()
    vec_n = ((vec2 - vec1) ** vec23)/((vec2 - vec1) ** vec23).norm()
    nbc = vec_n ** vec23
    buf = array([vec23[0], vec23[1], vec23[2], \
        nbc[0], nbc[1], nbc[2], \
        vec_n[0], vec_n[1], vec_n[2]])
    matrix = ndarray(shape=(3, 3), buffer=buf)
    return matrix


def build_coord(vec1, vec2, vec3, dist, angle, torsion, matrix=None):
    """"
    Builds coordinates assuming a reference frame:
    A (-1/-1/0)
    B (-1/0/0)
    C (0/0/0)
    """
    angle = math.radians(180 - angle)
    torsion = math.radians(torsion)
    # create initial vector
    vec_x = dist * math.cos(angle)
    vec_y = dist * math.cos(torsion) * math.sin(angle)
    vec_z = dist * math.sin(torsion) * math.sin(angle)

    # TODO: cache matrices for (dist/angle/torsion) (memoize pattern?)
    vec_d2 = Vector([vec_x, vec_y, vec_z])
    if matrix == None:
        matrix = get_ref_matrix(vec1, vec2, vec3)
    result_vec = vec_d2.right_multiply(matrix) + vec3
    return result_vec

