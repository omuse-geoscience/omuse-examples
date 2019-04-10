import cPickle
import os

def my_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--file", 
        dest="argfile", default = None,  
        help="arg file")
    return result

if __name__=="__main__":
    q=my_option_parser()
    o2,a2=q.parse_args()    
  
    try:
      label=a2[0] if a2 else o2.argfile
      if os.path.isfile(label):
        print label, "exists"
        print
      f=open(label,"r")
      kwargs=cPickle.load(f)
      for k in kwargs:
        print k,kwargs[k]
    except:
      raise
