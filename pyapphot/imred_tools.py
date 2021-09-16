''' All the reduction tasks are defined here which are called in image_management module.'''

from pyraf import iraf as ir
from iraf import imred
from iraf import ccdred
import numpy as np
from .utils import *
from astropy.io import fits

__all__ = ['master_bias_gen','master_dark_gen','master_flat_gen','normalize_frame','Pyccdproc','addpresuff2filename','winbin_fits','shift_gen']

def master_bias_gen(par, bias='@biaslist', combine='average', output='masterbias.fits',):
    if os.path.exists(output): os.remove(output)
    task = ir.zerocom
    task.combine = combine
    task.gain = par[0]
    task.rdnoise = par[1]
    task.ccdtype = ' '
    task(input=bias, output=output)
    return output

def master_dark_gen(par, dark='@darklist', combine='average', scale='none', output='masterdark.fits'):
    if os.path.exists(output): os.remove(output)
    task = ir.darkcom
    task.combine = combine
    task.gain = par[0]
    task.rdnoise = par[1]
    task.ccdtype = ''
    task.process = 'no'
    task.scale = scale
    task(input=dark, output=output)
    return output

def master_flat_gen(par, flat='@flatlist', combine='average', normalize=False, output='masterflat.fits'):
    # if os.path.exists(output): os.remove(folder+'/'+output)
    print(flat)
    task = ir.flatcom
    task.combine = combine
    task.gain = par[0]
    task.rdnoise = par[1]
    task.ccdtype = ' '
    task.process = 'no'
    task.subsets = 'no'
    task(input=flat, process='no', output=output)
    if normalize==True: output = normalize_frame(output, output)
    return output

def normalize_frame(frame, output):
    from astropy.io import fits
    f,h = fits.getdata(frame, header=True)
    m = np.max(f)
    fn = f/m
    if os.path.exists(output): os.remove(output)
    fits.writeto(output, fn, h)
    return output

def Pyccdproc(data='@datalist', bias='masterbias.fits', dark='', flat='masterflatn.fits', illum='', fringe='', output='bf_//@datalist', verbose=True):
    # ir.chdir(folder)
    for fn in filenames2list(output):
        if os.path.exists(fn): os.remove(fn)
    task = ir.ccdproc
    task.unlearn()
    task.ccdtype = ' '
    task.noproc = 'no'
    task.fixpix = 'no'
    task.overscan = 'no'
    task.trim = 'no'
    correclog = []
    if bias!='' and bias is not None:
        task.zerocor = 'yes'
        correclog.append('bias')
        task.zero = bias
    else: task.zerocor = 'no'
    if dark!='' and dark is not None:
        task.darkcor = 'yes'
        correclog.append('dark')
        task.dark = dark
    else: task.darkcor = 'no'
    if flat!='' and flat is not None:
        task.flatcor = 'yes'
        correclog.append('flat')
        task.flat = flat
    else: task.flatcor = 'no'
    if illum!='' and illum is not None:
        task.illumcor = 'yes'
        correclog.append('illum')
        task.illum = illum
    else: task.illumcor = 'no'
    if fringe!='' and fringe is not None:
        task.fringecor = 'yes'
        correclog.append('fringe')
        task.fringe = fringe
    else: task.fringecor = 'no'
    task.readcor = 'no'
    task.scancor = 'no'

    task.illum = illum
    task.fringe = fringe

    task(images=data, output=output)
    if verbose:
        correclog = ', '.join(correclog[:-1]) + ' and ' + correclog[-1]
        print(correclog + ' correction from data is done')

def addpresuff2filename(filename,pre='',suf=''):
    if type(filename)==str:
        d,b = os.path.dirname(filename),os.path.basename(filename)
        bn,be = os.path.splitext(b)
        return os.path.join(d,pre+bn+suf+be)
    if isinstance(filename,(list,tuple)):
        nfnames = []
        for fn in filename: nfnames.append(addpresuff2filename(fn,pre,suf))
        return nfnames

def winbin_fits(image,win,bin,bintype,output=''):
    if output!='': output = filenames2list(output)
    else: output = filenames2list(image)
    for i,fn in enumerate(filenames2list(image,forcelist=True)):
        d,h = fits.getdata(fn,header=True)
        if win!='full': d = d[win[2]-1:win[3]+win[2]-1,win[0]-1:win[0]+win[1]-1]
        dc = bintype(d.reshape(d.shape[0],-1,bin[0]),2)
        dr = bintype(dc.reshape(-1,bin[1],dc.shape[1]),1)
        if (type(output)==str and output==image) or output[i]==fn:fits.update(fn,dr,header=h)
        elif type(output)==str: fits.writeto(output,dr,header=h)
        elif isinstance(output,(list,tuple)): fits.writeto(output[i],dr,header=h)
    return output

def shift_gen(data, cooval, box):
    task = ir.imcntr
    task.cboxsize = box
    # print 'cooval:',
    # print cooval
    # print data
    task(input=data, x_init=cooval[0], y_init=cooval[1], Stdout='ctrfile')
    icoo = np.loadtxt('ctrfile', usecols=(2,4), unpack=False)
    return icoo[0] - icoo
