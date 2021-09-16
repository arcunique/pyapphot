from warnings import showwarning
import os
from . import joinpath2root as jpr

def j2s(*args):
    return jpr('scratch',*args)

def j2c(*args):
    return jpr('cache',*args)

class IrafWarning(Warning): pass
class PyapphotWarning(Warning): pass
class PyapphotException(Exception): pass

def warn(message,category=IrafWarning):
    showwarning(message,category,'','')

def filenames2list(filename, samepathstate=True, forcelist=False):
    import glob, os
    if filename is None: return None
    filenames = []
    if type(filename) == str and ',' not in filename and '*' not in filename and '@' not in filename: filenames = filename if not forcelist else [filename]
    elif type(filename)==str and '@' in filename:
        pref, suff = filename.split('@')
        while '//' in pref or '\\' in pref: pref = pref.replace('//','').replace('\\','')
        rawfile,suff = suff.split('//') if '//' in suff else (suff,'')
        for line in open(rawfile).readlines():
            filenames.append(pref+line.rstrip('\n')+suff)
    elif type(filename)==str and '*' in filename:
        if '[' in filename and filename.endswith(']'):
            filename,index = filename.split('[')
            index = '['+index
        else: index = ''
        for item in glob.glob(filename):
            filenames.append(item+index)
    elif not isinstance(filename, str) and hasattr(filename,'__len__'): filenames = filename
    elif type(filename)==str and ',' in filename: filenames = filename.split(',')
    if samepathstate: return sorted(filenames) if isinstance(filenames,(list,tuple)) else filenames
    if type(filenames)==str:
        if os.path.abspath(filenames)==filenames: return os.path.basename(filenames)
        return os.path.abspath(filenames)
    files = []
    for fn in filenames:
        if os.path.abspath(filenames)==fn: files.append(os.path.basename(fn))
        else: files.append(os.path.abspath(fn))
    return sorted(files)

def lendata(filename):
    if os.path.exists(filename):
        with open(filename) as fr: lines = fr.readlines()
        return sum(['#' not in l for l in lines])

def filterdict(dictobj:dict, refkey:str, axis:(int,None)=None ,values=None, lolim=None, hilim=None):
    import numpy as np
    newdict = {key: np.array(val) for key,val in dictobj.items() if key!='diminfo'}
    if 'diminfo' in dictobj: newdict['diminfo'] = dictobj['diminfo']
    if axis is None:
        if 'diminfo' not in dictobj: axis=0
        else: axis = dictobj['diminfo'][refkey]
    refval = newdict[refkey].copy()
    condition = True
    if values is not None:
        condition &= np.in1d(refval,values)
    if lolim is not None:
        condition &= refval>=lolim
    if hilim is not None:
        condition &= refval<=hilim
    if not np.all(condition):
        for key,val in newdict.items():
            if key=='diminfo': continue
            if np.ndim(val)==1 and len(val)==len(refval):
                if 'diminfo' not in dictobj or newdict['diminfo'][key]==axis:
                    newdict[key] = val[condition]
            elif np.ndim(val)!=1 and np.ndim(val)>=axis:
                if axis==0: newdict[key] = val[condition]
                else:
                    shape = list(val.shape)
                    ind = [np.ones(s, dtype=bool) for s in shape]
                    ind[axis] = condition
                    shape[axis] = sum(condition)
                    index = True
                    for i in ind: index = np.tensordot(index,i,axes=0)
                    newdict[key] = val[index].reshape(shape)
    return newdict

def store_filenames(files, reject=(), saveto=''):
    files = filenames2list(files, forcelist=True)
    with open(saveto, 'w') as fw:
        fw.writelines([file+'\n' for file in files if file not in reject])

def segregate(x, shape):
    i = 0
    y = []
    for sh in shape:
        y.append(x[i:i+sh])
        i = i+sh
    return y