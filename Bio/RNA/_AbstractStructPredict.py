class AbstractStructPredict(object):
    """An abstract class for RNA secondary structure prediction.
    All inheriting classes should implement the predict() method that
    returns a Secstruc class
    """
    
    def __init__(self):
        """
        """
        pass
    def predict(self, RNA):
        return None
