"""Tools for making calculation results persistent.
"""
import os,pickle

class persistent(object):
    """Decorator for making the results persistent.
    """

    def __init__(self,function):
        self.function=function

    def __call__(self,*args,**kwargs):
        if "persistent_file" in kwargs:
            persistent_file=kwargs.pop("persistent_file")
            if os.path.exists(persistent_file): #read from the file
                result = pickle.load(open(persistent_file,"r"))
            else: #create it, write it
                result = self.function(*args,**kwargs)
                thefile = open(persistent_file,"w")
                pickle.dump(result,thefile)
                thefile.close()
            return result
        else:
            return self.function(*args,**kwargs)

