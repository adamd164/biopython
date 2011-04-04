#!/usr/bin/env python
#
# test_coord_builder.py
#
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""
Tests for the CoordBuilder module.
"""

__author__ = "Kristian Rother, Magdalena Musielak, Tomasz Puton"
__license__ = "Biopython License"
__credits__ = ["Janusz Bujnicki"]
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from unittest import main, TestCase
from Bio.PDB.CoordBuilder import build_coord
from Bio.PDB.Vector import Vector, calc_angle, calc_dihedral
import math
from random import random

class CoordBuilderTests(TestCase):
    """
    Tests for the CoordBuilder functions.
    """
    def setUp(self):
        """Loads the A residue to start with."""
        # coordinates for A,B,C in the reference frame
        self.ref_coord = [
            Vector([-1.0, -1.0, 0.0]),
            Vector([-1.0,  0.0, 0.0]),
            Vector([ 0.0,  0.0, 0.0])
            ]
        # some other coordinates for A,B,C
        self.sample_coord = [
            Vector([1.0, 0.0, 0.0]),
            Vector([0.0, 0.0, 0.0]),
            Vector([0.0, 1.0, 0.0])
            ]

    def _get_random(self, dist=None, angle=None, torsion=None):
        """Returns random distance, angle, and torsion"""
        if dist == None:
            dist = random()*5.0
        if angle == None:
            angle = random()*180
        if torsion == None:
            torsion = 180-random()*360.0
        return dist, angle, torsion
            
    def _measure(self, vec1, vec2, vec3, vec4):
        """Returns a dist,angle,torsion tuple for the given four coords."""
        dist = (vec3 - vec4).norm()
        angle = math.degrees(calc_angle(vec2, vec3, vec4))
        torsion = math.degrees(calc_dihedral(vec1, vec2, vec3, vec4))
        return dist, angle, torsion

    def _test_coord_samples(self, coord, examples):
        """
        Runs the coord building on a set of
        dist, angle, torsion values.
        """
        for dis, ang, tor in examples:
            vec1, vec2, vec3 = coord
            vec4 = build_coord(vec1, vec2, vec3, dis, ang, tor)
            result = self._measure(vec1, vec2, vec3, vec4)
            self.assertDistance(result, dis)
            self.assertAngle(result, ang)
            self.assertTorsion(result, tor)

    def assertDistance(self, got, expected):
        """Assertion of a distance value."""
        self.assertAlmostEqual(got[0], expected, 3)

    def assertAngle(self, got, expected):
        """Assertion of an angle."""
        self.assertAlmostEqual(got[1], expected, 3)

    def assertTorsion(self, got, expected):
        """Assertion of a torsion angle."""
        self.assertAlmostEqual(got[2], expected, 3)

        
    def test_build_coord(self):
        """Checks whether the P+O5' atoms are constructed."""
        vec1, vec2, vec3 = self.sample_coord
        vec4 = build_coord(vec1, vec2, vec3, 1.0, 90.0, 0)
        self.assertTrue(isinstance(vec4, Vector))

    def test_distance(self):
        """Checks whether the distance is OK."""
        vec1, vec2, vec3 = self.sample_coord
        vec4 = build_coord(vec1, vec2, vec3, 2.34, 90.0, 0.0)
        result = self._measure(vec1, vec2, vec3, vec4)
        self.assertDistance(result, 2.34)

    def test_angle(self):
        """Checks whether the angle is OK."""
        vec1, vec2, vec3 = self.sample_coord
        vec4 = build_coord(vec1, vec2, vec3, 1.0, 90.0, 0.0)
        result = self._measure(vec1, vec2, vec3, vec4)
        self.assertAngle(result, 90.0)

    def test_dihedral(self):
        """Checks whether the torsion angle is OK."""
        vec1, vec2, vec3 = self.sample_coord
        vec4 = build_coord(vec1, vec2, vec3, 1.0, 90.0, 90.0)
        result = self._measure(vec1, vec2, vec3, vec4)
        self.assertTorsion(result, 90.0)

    def test_multiple(self):
        """Runs the procedure on known examples."""
        self._test_coord_samples(self.sample_coord, EXAMPLES)
        self._test_coord_samples(self.ref_coord, EXAMPLES)

    def test_dist_power(self):
        """Checks the distance with random values."""
        random_examples = [self._get_random(angle=0.0, torsion=0.0) \
            for x in range(100)]
        self._test_coord_samples(self.sample_coord, random_examples)

    def test_angle_power(self):
        """Checks the angle with random values."""
        random_examples = [self._get_random(dist=1.51, torsion=0.0) \
            for x in range(100)]
        self._test_coord_samples(self.sample_coord, random_examples)

    def test_torsion_power(self):
        """Checks the torsion with random values."""
        random_examples = [self._get_random(dist=1.51, angle=32.0) \
            for x in range(100)]
        self._test_coord_samples(self.sample_coord, random_examples)
            
    def test_max_power(self):
        """Runs the procedure many times with random values."""
        random_examples = [self._get_random() for x in range(100)]
        self._test_coord_samples(self.sample_coord, random_examples)


EXAMPLES = [
    ( 2.0, 90.0, 0.0),
    ( 1.0, 45.0, 0.0),
    ( 1.0, 45.0, -90.0),
    ( 1.0, 45.0, 90.0),
    ( 1.0, 45.0, 180.0),
    ( 1.0, 45.0, 135.0),
    ]


if __name__ == '__main__':
    main()
  
