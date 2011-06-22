#Changed by Asaf Peer 2011
import RNAalifoldPredictor


_FormatToPredictors={'RNAalifold' : RNAalifoldPredictor.AFPredictor}
def predictStruct(seq,format):
    if format not in _FormatToPredictors:
        raise ValueError("choose one of %s"%str(_FormatToPredictors.keys()))
    pred_class = _FormatToPredictors[format]
    return pred_class().predict(seq)

    
    

    
        
