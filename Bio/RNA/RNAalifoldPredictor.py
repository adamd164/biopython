#Copyright Asaf Peer 2011
import tempfile
from Bio.RNA._AbstractStructPredict import AbstractStructPredict
from Bio import AlignIO
from Bio.Vienna.Applications import RNAalifoldCMD
from Secstruc import Secstruc
class AFPredictor(AbstractStructPredict):
    """Predict the secondary structure of a multiple sequence alignment using RNAalifold.
    Use the Vienna module.
    """
    
    def __init__(self):
        """
        """
        pass

    def predict(self, MSA):
        """Predict the consensus structure of a MSA.
        
        Arguments:
        - `MSA`: A MultipleSequenceAlignment object
        """
        fileh = tempfile.NamedTemporaryFile()
        AlignIO.write(MSA,fileh,'clustal')
        fileh.flush()
        alifold = RNAalifoldCMD(filename=fileh.name,noPostScript=True)
        stdout, stderr=alifold()
        #Assume RNAalifold return 2 lines, the first contains the consensus sequence
        #and the second the consensus fold
        lines=stdout.split('\n')
        import sys
        sys.stderr.write(str(lines))
        fold = lines[1].split()[0]
        return Secstruc(fold)
    
        

        
