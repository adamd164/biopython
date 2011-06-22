#Copyright 2011 Asaf Peer
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for interacting with the Vienna RNA package
These classes follow the AbstractCommandline interfaces for running
programs.
"""

from Bio.Application import _Option, _Switch, AbstractCommandline

class _ViennaMinimalCommandLine(AbstractCommandline):
    """Base Commandline for Vienna package

    This class deals with shared options:

     - C             Constraints file
     - T             Temperature
     - d             dangling ends (0|1|2|3)
     - noLP          without lonely-pairs
     - noGU          Do not allow GU pairs
     - nocloseGU     Do not allow GU pairs at the end of helices
     - P             Read parameters from paramfile
     - nsp           Allow other pairs
     - S             Scale the energy
     - circ          Assume circular RNA
     - noPS          Don't produce postscript drawing of the structure
     
    """
    
    def __init__(self, cmd=None, **kwargs):
        extra_parameters = [
            _Option(['-C','constrains'],'''Calculate structures subject to constraints.   The  program  reads  first  the
              sequence,  then  a string containing constraints on the structure encoded with
              the symbols: | (the corresponding base  has  to  be  paired  x  (the  base  is
              unpaired) < (base i is paired with a base j>i) > (base i is paired with a base
              j<i) and matching brackets ( ) (base i pairs base j)  With  the  exception  of
              "|", constraints will disallow all pairs conflicting with the constraint. This
              is usually sufficient to enforce the constraint, but occasionally a  base  may
              stay  unpaired in spite of constraints. PF folding ignores constraints of type
              "|".
        ''',
                    filename=True,equate=False),
            _Option(['-T','temperature'],'Rescale energy parameters to a temperature of temp C. Default is 37C.',
                    equate=False),
            _Option(['-P','paramfile'],'Read energy parameters from paramfile, instead of using the default  parameter set.',
                    equate=False,filename=True),
            _Option(['-S','scale'],'In the calculation of the pf use scale*mfe as an  estimate  for  the  ensemble free  energy (used to avoid overflows)',
                    equate=False),
            _Option(['-nsp','pairs'],'Allow other pairs',equate=False),
            _Switch(['-d0','dang-ignore'],'Ignore dangling ends'),
            _Switch(['-d1','dang-unpaired'],'only unpaired bases  can  participate  in  at most one dangling end, this is the default for mfe folding but unsupported for the partition function folding'),
            _Switch(['-d2','dang-helix'],'dangling energies  will  be  added  for  the bases adjacent to a helix on both sides in any case'),
            _Switch(['-d3','dang-coaxial'],'mfe  folding  will  allow  coaxial  stacking of adjacent helices in multi-loops'),
            _Switch(['noLP','noLonelyPairs'],'Produce structures without lonely pairs (helices of length 1)'),
            _Switch(['-noGU','noGU'],'Do not allow GU pairs'),
            _Switch(['-noCloseGU','noCloseGU'],'Do not allow GU pairs at the end of helices'),
            _Switch(['-circ','circular'],'Assume  a  circular (instead of linear) RNA molecule'),
            _Switch(['-noPS','noPostScript'],'Do not produce postscript drawing of the mfe structure')
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)
                    
class RNAalifoldCMD(_ViennaMinimalCommandLine):
    """
    Run RNAalifold on multiple sequence alignment.
    """
    def __init__(self, cmd="RNAalifold", **kwargs):
        self.parameters = [
            _Option(['','filename'],'The MSA file',equate=False,filename=True),
            _Option(['-cv','covariance'],'Set  the  weight  of  the  covariance  term  in the energy function to factor',
                    equate=False),
            _Option(['-nc','non-compatible'],'Set the penalty for non-compatible sequences in the  covariance  term  of  the energy function to factor',
                    equate=False),
            _Switch(['-mis','mis'],'Output  "most informative sequence" instead of simple consensus'),
            _Switch(['-E','endgaps'],'Score pairs with endgaps same as gap-gap pairs'),
            _Switch(['-p','partition'],'Calculate  the partition function and base pairing probability matrix'),
            _Switch(['-color','color'],'Produce a colored version of the consensus strcture plot "alirna.ps"'),
            _Switch(['-aln','alignment'],'Produce  a  colored  and structure annotated alignment in PostScript format in the file "aln.ps" in the current directory')
            ]

        _ViennaMinimalCommandLine.__init__(self, cmd, **kwargs)
        
