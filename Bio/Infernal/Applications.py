#Copyright 2011 Aasf Peer.
"""
Code for executing Infernal
"""
#Only options and switches that appear in the tutorial are listed here
#Should add more (if I'll need it sometimes...)
from Bio.Application import _Option,_Switch,_Argument,AbstractCommandline

class cmAlignCMD(AbstractCommandline):
    def __init__(self, cmd="cmalign", **kwargs):
        self.parameters = [
            _Option(['-o','outfile'],'Print the MSA to a file',equate=False),
            _Option(['--withali','withali'],'Include the alignment used to '
                    'generate the model (if --rf of --gapthresh were used to'
                    'generate the model they should be used here as well',equate=False),
            _Option(['--gapthresh','gapthresh'],'Threshold for percent of gaps'
                    'in column to be considered as insertion',
                    equate=False),
            _Switch(['--rf','rf'],'Uses the RF line to build a model '
                    '(useful when expanding seed alignment that has a GF line from cmalign)'),
            _Switch(['-l','local'],'Performs a local alignment'),
            _Switch(['--sub','sub'],'Some sequences are not complete'),
            _Switch(['-p','posterior'],'Print the posterior probability of the alignment'),
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            _Argument(['fastafile'],'fasta file',filename=True,is_required=True)

            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class cmBuildCMD(AbstractCommandline):
    def __init__(self, cmd="cmbuild", **kwargs):
        self.parameters = [

            _Option(['--gapthresh','gapthresh'],'Threshold for percent of gaps'
                    'in column to be considered as insertion',
                    equate=False),
            _Option(['--cmaxid','cmaxid'],'Maxmimum identity between sequences'
                    'in different clusters',equate=False),
            _Option(['--ctarget','ctarget'],'Number of clusters',equate=False),
            _Option(['--cdump','cdump'],'Dump the clusters into a file',equate=False),
            _Switch(['--corig','corig'],'Create a model from all the sequences'
                    'in addition to the clustered models'),
            _Switch(['--call','call'],'Create a model from each sequence in the sto file'),
            _Option(['--refine','refine'],'Refine the model by iteratively '
                    'generating a model and aligning the sequences to it, the '
                    'alignment is saved to the given filename',equate=False),
            _Switch(['--gibbs','gibbs'],'Gibbs sample from the refined models'
                    '(only with --refine)'),
            _Switch(['--rf','rf'],'Uses the RF line to build a model '
                    '(useful when expanding seed alignment that has a GF line from cmalign)'),
            _Option(['--rsearch','rsearch'],'build RSEARCH model using the '
                    'given RIBOSUM matrix (use only with --call or a single'
                    'input sequence',equate=False),
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            _Argument(['stkfile'],'multiple sequence alignment in Stockholm format',
                    filename=True,is_required=True)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class cmCalibrateCMD(AbstractCommandline):
    def __init__(self, cmd="cmcalibrate", **kwargs):
        self.parameters = [

            _Option(['--forecast','forecast'],'Forecast the time of execution, not calibrating',
                    equate=False,is_required=False),
            _Argument(['modelname'],'name of model',filename=True,is_required=True)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class cmEmitCMD(AbstractCommandline):
    #Options not impolemented
    def __init__(self, cmd="cmemit", **kwargs):
        self.parameters = [
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            _Argument(['outputfile'],'sequence output file',filename=True,is_required=True)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class cmScoreCMD(AbstractCommandline):
    #Options not impolemented
    def __init__(self, cmd="cmscore", **kwargs):
        self.parameters = [
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        
class cmSearchCMD(AbstractCommandline):
    def __init__(self, cmd="cmsearch", **kwargs):
        self.parameters = [
            _Option(['--forecast','forecast'],'Forecast the time of execution, not searching',
                    equate=False,is_required=False),
            _Option(['-Z','size'],'Database size in MB',equate=False),
            _Switch(['-g','glocal'],'Performs a glocal alignment search',equate=False),
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            _Argument(['database'],'Database file', filename=True,is_required=True)

            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

class cmStatCMD(AbstractCommandline):
    #Options not implemented
    def __init__(self, cmd="cmstat", **kwargs):
        self.parameters = [
            _Argument(['modelname'],'name of model',filename=True,is_required=True),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)



    
