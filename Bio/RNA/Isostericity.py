#!/usr/bin/env python
#
# Isostericity.py 
#
"""
Module managing RNA isostericity matrices.

TODO: add reference
"""
__author__  = "Pawel Skiba"
__contributors__ = "Kristian Rother, Magdalena Musielak"
__credits__ = ["Tomasz Puton", "Janusz M. Bujnicki"]
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Production"

import os

DATA_PATH = os.path.split(__file__)[0]
ISO_MATRIX_DATA = DATA_PATH + os.sep + "isostericity_matrices.txt"
ISO_DISCREPANCY_CUTOFF = 1.0

class IsostericityError(Exception): 
    pass

class IsostericityMatrices:
    """
    Isostericity matrices implementation in Python.
    
    In the isostericity matrices, base pairs are represented by
    two-character strings 'AU', 'GC', 'AA', etc.
    Base pair types are represented by the abbreviated Leontis-Westhof
    nomenclature 'cWW', 'tWH', 'cSS' etc.
    """
    def __init__(self, data_file=ISO_MATRIX_DATA):
        """Reads isostericity matrices from a data file."""
        self.matrices = {}
        self.parse_matrices(open(data_file))
        
    def parse_matrices(self, handle):
        """Parses a file with isostericity matrices into a dictionary."""
        pair1 = None
        bp_type = None
        for line in handle:
            if line.startswith("Type: "):
                data = line.split(" ")
                pair1 = data[1]
                bp_type = data[2].strip()
                self.matrices.setdefault(pair1, {})
                self.matrices[pair1].setdefault(bp_type, {})
            elif pair1 and not line.startswith("\n"):
                data = line.split(": ")
                pair2 = data[0]
                if data[1] == '\n':
                    iso_discrepancy = None
                else:
                    iso_discrepancy = float(data[1].strip()) 
                self.matrices[pair1][bp_type][pair2] = iso_discrepancy
        
    def get_isostericity(self, bp1, bp2, bp_type):
        """
        Returns the isostericity value for replacing basepair1 by basepair2
        with basepair1 having the given type.
        """
        try:
            result = self.matrices[bp1][bp_type][bp2] 
        except KeyError:
            raise IsostericityError("No value in isostericity matrices for: "\
                +bp1+"->"+bp2+" ("+bp_type+")")
        return result
        
    def check_isostericity(self, bp1, bp2, bp_type, \
                           max_discrepancy=ISO_DISCREPANCY_CUTOFF):
        """
        Returns True if basepair1 can be replaced by an isosteric basepair2 
        with basepair1 having the given type, and
        the discrepancy is smaller than the given value
        interaction type is interact_type
        """
        result = self.get_isostericity(bp1, bp2, bp_type)
        return result <= max_discrepancy and result != None
        
    def show_isosteric_bp(self, bp1, bp_type, \
                          max_discrepancy=ISO_DISCREPANCY_CUTOFF):
        """
        Returns all base pairs isosteric to bp1 with the given base pair type
        and maximum discrepancy.
        """
        pairs = self.matrices[bp1][bp_type]
        result = []
        for basepair in pairs:
            if pairs[basepair] <= max_discrepancy and pairs[basepair] != None:
                result.append(basepair)
        return tuple(result)
        
