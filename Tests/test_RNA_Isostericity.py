#!/usr/bin/env python
#
# test_RNA_Isostericity.py 
#
"""
Unit Tests for the module managing RNA isostericity matrices.
"""

__author__  = "Pawel Skiba"
__contributors__ = "Kristian Rother, Magdalena Musielak"
__credits__ = ["Tomasz Puton", "Janusz M. Bujnicki"]
__maintainer__ = "Pawel Skiba"
__email__ = "pw.skiba@gmail.com"
__status__ = "Production"


from unittest import main, TestCase
from Bio.RNA.Isostericity import IsostericityMatrices

class IsostericityMatricesTests(TestCase):
    """Tests for the module parsing isostericity matrices."""
    
    def test_check_isostericity(self):
        """Isostericity has a default cutoff."""
        ima = IsostericityMatrices()
        result = ima.check_isostericity('AC', 'GU', 'cWW')
        self.assertTrue(result)
        result2 = ima.check_isostericity('AA', 'GG', 'cWW')
        self.assertFalse(result2)
        result3 = ima.check_isostericity('AC', 'UU', 'cWW')
        self.assertFalse(result3)

    def test_check_isostericity_cutoff(self):
        """Setting the cutoff influences isostericity results."""
        ima = IsostericityMatrices()
        result = ima.check_isostericity('UU', 'AA', 'tWH',  2.0)
        self.assertFalse(result)
        result2 = ima.check_isostericity('UU', 'AA', 'tWH',  3.0)
        self.assertTrue(result2)
        
    def test_check_isostericity_none(self):
        """Some isostericity values may be empty."""
        ima = IsostericityMatrices()
        result = ima.check_isostericity('AA', 'GA', 'tWS')
        self.assertFalse(result)

    def test_get_isosteric_pairs(self):
        """A list of isosteric pairs can be requested."""
        ima = IsostericityMatrices()
        result = ima.get_isosteric_pairs('AC', 'cWW')
        self.assertEqual(result, ('AC', 'GU'))
   
    def test_get_isosteric_pairs_cutoff(self):
        """Isosteric pairs can be influenced by cutoff."""
        ima = IsostericityMatrices()
        result = ima.get_isosteric_pairs('AC', 'cWW', 3.0)
        result = list(result)
        result.sort()
        self.assertEqual(result, ['AC', 'AU', 'CG', 'GC', 'GU', 'UA', 'UU'])
   
if __name__ == '__main__':
    main()
    
