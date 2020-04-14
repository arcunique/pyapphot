import os
rootdir = __path__[0]
from .version import __version__

def joinpath2root(*args):
    if len(args)==1 and os.path.exists(args[0]): return args[0]
    return os.path.join(rootdir,*args)

configpath = joinpath2root('config')