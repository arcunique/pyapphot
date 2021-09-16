''' This module allows examination of the images, reduction of the raw frames and alignment of frames. '''

import logging
import numpy as np
from . import imexam_mod as imexam
import matplotlib
from .utils import *
from .imred_tools import *
from . import _configpath
from astropy.io import fits
from configparser import ConfigParser
# matplotlib.use('Qt5Agg')
imexam.set_logging(logging.CRITICAL)
__all__ = ['get_psfcentroid','imexamine','imcorrection','alignimages','save_matched_coords']

try:
    from pyraf import iraf as ir
    if not os.path.exists('login.cl'): os.system('mkiraf -t xgterm -i > /dev/null')
    irafimported = True
except: irafimported = False


def get_psfcentroid(image,x,y,box=20):
    '''
    calculates psfcentroid using the iRAF imcntr task

    :param image: Image name or list of images: str or list
    :param x: x co-ordinates: int or float
    :param y: y co-ordinate: int or float
    :param box: box size, within which it will find the centroid: int or float
    :return: centroid (x,y): numpy.ndarray
    '''
    if not irafimported:
        warn('pyraf/IRAF not installed. Centering algorithm will not work.')
        return
    image = filenames2list(image)
    if type(image)==str:
        task = ir.imcntr
        task.cbox = box
        if os.path.exists(j2s('__coo__')): os.remove(j2s('__coo__'))
        task(input=image, x_init=x, y_init=y, Stdout=j2s('__coo__'))
        res = np.loadtxt(j2s('__coo__'), usecols=(2,4), unpack=False)
        return res
    coo = []
    for im in image: coo.append(get_psfcentroid(im,x,y,box))
    return np.array(coo)

def gaussfit_psf(image,x,y,data=None,delta=50,amp2bg_cutoff=0.5,):
    '''
    Performs 2D Gaussian fit on the PSF around x,y and return the xcenter,ycenter and uncertainties.
    :param image: Image name: str
    :param x: x co-ordinates: int or float
    :param y: y co-ordinate: int or float
    :param data: Array of data stored in the sci extension of the fits file around x,y : numpy.ndarray
    :param delta: box size, within which it will find the centroid: int or float
    :param amp2bg_cutoff: amplitude to background ratio which will identify whether it is part of the PSF or background sky: float/int
    :return: xcenter, ycenter, uncertainty in xcenter and uncertainty in ycenter
    '''
    from imexam.imexamine import Imexamine as imexa
    # imexam.set_logging(logging.CRITICAL)
    image = filenames2list(image)
    if type(image)==str:
        if not data:
            if '[' in image and image.endswith(']'):
                imname, sci = image.split('[')[0], int(image.split('[')[1].rstrip(']'))
            else: imname, sci = image, 0
            data = fits.getdata(imname,sci)
        amp, x, y, xsig, ysig = imexa().gauss_center(x, y, data=data, delta=delta)
        delta = int(delta)
        xx = int(x)
        yy = int(y)
        background = data[yy - delta:yy + delta, xx - delta:xx + delta].min()
        if amp / background > amp2bg_cutoff:
            return amp, x, y, xsig, ysig

def get_coords_filename_from_prefix(image, prefix=''):
    image = filenames2list(image, forcelist=False)
    if type(image)==str:
        dirn,filen = os.path.dirname(image), os.path.basename(image)
        filen,fileext = os.path.splitext(filen)
        return os.path.join(dirn,prefix+filen+'.coo')
    coofs = []
    for im in image:
        coofs.append(get_coords_filename_from_prefix(prefix=prefix,image=im))
    return coofs

def get_image_assigned2coofile(coofile,image=''):
    if not os.path.exists(coofile) or os.stat(coofile).st_size==0: return
    lines = open(coofile).readlines()
    if '#' in lines[0]:
        if image=='': return lines[0].split()[1]
        elif '.fits' in lines[0].split()[1]: return lines[0].split()[1]==image

def set_image_assigned2coofiles(coofile,image):
    if not os.path.exists(coofile) or os.stat(coofile).st_size==0: return
    lines = open(coofile).readlines()
    if not get_image_assigned2coofile(coofile) and lendata(coofile):
        lines.insert(0,'# '+image+'\n')
        with open(coofile,'w') as fw: fw.writelines(lines)
        return
    if not get_image_assigned2coofile(coofile,image):
        lines[0] = '# '+image+'\n'
        with open(coofile,'w') as fw: fw.writelines(lines)


class imexamine(object):

    def __init__(self, set_xpa_method_local=False):
        '''
        contains methods to examine the images.
        :param set_xpa_method_local: This is a flag which sets the method to be 'local' if True else 'inet': bool
        '''
        self.viewer = None
        self.xpa = None
        self._localxpa = False
        if set_xpa_method_local:
            self._org_xpa_method = os.getenv('XPA_METHOD')
            if self._org_xpa_method!='local':
                os.environ['XPA_METHOD'] = 'local'
            self._localxpa = True
        else:
            self._org_xpa_method = os.getenv('XPA_METHOD')
            if self._org_xpa_method != 'inet':
                os.environ['XPA_METHOD'] = 'inet'

    def openviewer(self,ds9_target=None,wait_time=10):
        '''
        opens the ds9 window
        :param ds9_target: If ds9 is opened already, the title can be passed here to connect to the same ds9 window.
        :param wait_time: This is a value which keeps the processes waiting until the ds9 window is open, otherwise
        throws error as the process cannot find the ds9 target: int
        :return: None
        '''
        self.viewer = imexam.connect(target=ds9_target,viewer='ds9',wait_time=wait_time)
        # self.viewer.exam.setlog(level=logging.CRITICAL)

    def display(self,image, **kwargs):
        '''
        Displays the image
        :param image: image name: str
        :param kwargs: kwargs in openviewer
        :return: None
        '''
        if self.viewer is None: self.openviewer(**kwargs)
        self.viewer.load_fits(image)
        self.viewer.scale()
        self.viewer.zoomtofit()
        self.image_ondisplay = image
        self._already_displayed = True
        if self._localxpa: self.viewer.window.set_iraf_display()

    @property
    def next_image(self):
        self._imdeque.rotate(-1)
        if self._coodeque: self._coodeque.rotate(-1)
        return self._imdeque[0]

    @property
    def previous_image(self):
        self._imdeque.rotate(1)
        if self._coodeque: self._coodeque.rotate(1)
        return self._imdeque[0]

    @property
    def current_image(self):
        return self._imdeque[0]

    def mark_coords(self,image=None,coofile='',**kwargs):
        if image is not None: self.display(image, **kwargs)
        if coofile == '' or not os.path.exists(coofile) or os.stat(coofile).st_size == 0: return
        dat = np.loadtxt(coofile, ndmin=2)
        self.viewer.mark_region_from_array(np.c_[dat, np.arange(1, dat.shape[0] + 1)], textoff=90, size=60)
        return np.c_[dat, np.arange(1, dat.shape[0] + 1)]

    @property
    def _get_coofn_from_prefix(self):
        return get_coords_filename_from_prefix(image=self.image_ondisplay, prefix=self._coordsprefix)

    @property
    def get_cur_coofile(self):
        if self._singlecoo:
            if not self._im2coo: self._im2coo = get_image_assigned2coofile(self._singlecoo)
            else: set_image_assigned2coofiles(self._singlecoo, self._im2coo)
            if not self._im2coo or self._im2coo==self.image_ondisplay: return self._singlecoo
        if self._coodeque:
            cooass = self._coodeque[0]
            cooassim = get_image_assigned2coofile(cooass)
            if cooassim and cooassim != self.image_ondisplay: os.remove(cooass)
            return cooass
        if self._coordsprefix is not None:
            return self._get_coofn_from_prefix

    def clear_coords_file(self,coordsprefix='',image=None,):
        import glob
        if type(image)==str and image=='__all__':
            for fn in glob.glob(coordsprefix+'*.coo'):
                if os.path.exists(fn): os.remove(fn)
        else:
            coofile = self._get_coo_fn(prefx=coordsprefix,image=image)
            if type(coofile)==str and os.path.exists(coofile): os.remove(coofile)
            elif type(coofile)==list:
                for cf in coofile:
                    if os.path.exists(cf): os.remove(cf)

    def get_psfcenter(self,x,y,delta=50,amp2bg_cutoff=0.5,full_output=False):
        imdata = self.viewer.window.get_data().astype(float)
        # print(imdata)
        amp, x, y, xsig, ysig = self.viewer.exam.gauss_center(x, y, data=imdata, delta=delta)
        delta = int(delta)
        xx = int(x)
        yy = int(y)
        background = imdata[yy - delta:yy + delta, xx - delta:xx + delta].min()
        # print(amp,background)
        if amp/background>amp2bg_cutoff:
            if not full_output: return x, y
            return amp, x, y, xsig, ysig

    def examine_frames(self, image=None, coordfile='', markcoords=True, verbose=False, clear_prev_coords=False, **kwargs):
        '''
        Examines the frames depending on keys inserted.
        :param image: image name(s): str, str with '*' or '@', list, tuple
        :param coordfile (optional): File(s) where selected co-ordinates are to be stored. If this is skipped it will look for the 'coordprefix'
        key word which gets added before the filename and replaces the .fits extension with .coo extension. By default coordsprefix is '': str,
         str with '*' or '@', list, tuple
        :param markcoords: If set it marks the selected co-ordinates: bool
        :param verbose: bool
        :param clear_prev_coords: bool
        :param kwargs:
        :return:
        '''
        from collections import deque
        import matplotlib.pyplot as plt
        if image is None and self.viewer is None: raise PyapphotException('No ds9 window opened. '
       'At least open a ds9 window and also to run this command successfully either load an image prior or pass the images as an argument for \'images\'.')
        if image is not None:
            ds9_target = kwargs.pop('ds9_target',None)
            wait_time = kwargs.pop('wait_time',10)
            self._imdeque = deque(filenames2list(image,forcelist=True))
            self.display(self.current_image,ds9_target=ds9_target,wait_time=wait_time)
            self._already_displayed = False
        self._coodeque = None
        self._singlecoo = None
        if coordfile:
            coordfile = filenames2list(coordfile)
            if type(coordfile) == str:
                if clear_prev_coords:
                    if os.path.exists(coordfile): os.remove(coordfile)
                self._singlecoo = coordfile
                self._im2coo = kwargs.pop('imageassigned2coo',None)
            else:
                self._coodeque = deque(coordfile)
                if clear_prev_coords:
                    for coof in coordfile:
                        if os.path.exists(coof): os.remove(coof)
        else:
            self._coordsprefix = kwargs.pop('coordsprefix',None)
            if clear_prev_coords and self._coordsprefix:
                for coof in get_coords_filename_from_prefix(image,self._coordsprefix):
                    if os.path.exists(coof): os.remove(coof)
        cooass = self.get_cur_coofile
        if markcoords and cooass: self.mark_coords(coofile=cooass)
        while True:
            x, y, key = self.viewer.readcursor()

            if key == 'i':
                if irafimported:
                    print('entering iraf-imexamine mode, press q to exit')
                    ir.imexa()
                    print('exiting iraf-imexamine mode')
                else:
                    warn('pyraf/IRAF not installed. Cannot enter iraf-imexamine mode.')

            elif key == 's':
                cooc = get_psfcentroid(self.image_ondisplay,x,y)
                cooc = self.get_psfcenter(*cooc)
                if not cooass and self._singlecoo: print(self.image_ondisplay+' not assigned to '+coordfile)
                elif not cooass: print('no coordfile assigned')
                if np.all(cooc) and cooass:
                    x,y = cooc
                    if not os.path.exists(cooass) or not np.any(['#' not in l for l in open(cooass).readlines()]):
                        np.savetxt(cooass, [[x, y]], fmt='%0.02f')
                    else:
                        data = np.loadtxt(cooass, ndmin=2)
                        for coo in data:
                            if (coo[0] - x) ** 2 + (coo[1] - y) ** 2 < 400: break
                        else: np.savetxt(cooass, np.c_[data.T, [x, y]].T, fmt='%0.02f')
                    if coordfile:
                        set_image_assigned2coofiles(cooass,self.image_ondisplay)
                        print(self.image_ondisplay+' assigned to '+cooass)
                    if markcoords: self.mark_coords(image=self.image_ondisplay, coofile=cooass)

            elif key=='comma':
                cooc = get_psfcentroid(self.image_ondisplay, x, y)
                res = self.get_psfcenter(*cooc,full_output=True)
                if res: print(f'XC={res[1]}\tYC={res[2]}\nXC_err={res[3]}\tYC_err={res[4]}')
                else: print('No object found')

            elif key=='d':
                cooc = get_psfcentroid(self.image_ondisplay, x, y)
                cooc = self.get_psfcenter(*cooc)
                if np.all(cooc) and cooass:
                    x, y = cooc
                    if os.path.exists(cooass) and np.any(['#' not in l for l in open(cooass).readlines()]):
                        data = np.loadtxt(cooass, ndmin=2)
                        delind = []
                        for i in range(len(data)):
                            # print(data[i],x,y)
                            if (data[i][0]-x)**2+(data[i][1]-y)**2<400: delind.append(i)
                        # print(delind)
                        np.savetxt(cooass,np.delete(data,delind,axis=0),fmt='%0.02f')
                        set_image_assigned2coofiles(cooass, self.image_ondisplay)
                        if markcoords:
                            self.mark_coords(image=self.image_ondisplay,coofile=cooass)

            elif key == 'r':
                cooc = get_psfcentroid(self.image_ondisplay,x,y)
                if cooc is not None:
                    if not verbose: self.viewer.exam.log.setLevel(logging.CRITICAL)
                    # self.viewer.set_plot_pars('r','func','Moffat1D')
                    self.viewer.exam.radial_profile(cooc[0], cooc[1], data=self.viewer.window.get_data().astype(float), )
                    plt.pause(0.5)

            elif key=='n':
                self.display(self.next_image, ds9_target=ds9_target, wait_time=wait_time)
                cooass = self.get_cur_coofile
                if markcoords and cooass: self.mark_coords(coofile=cooass)

            elif key=='N':
                self._imdeque = deque(filenames2list(image, forcelist=True))
                self.display(self.previous_image, ds9_target=ds9_target, wait_time=wait_time)
                cooass = self.get_cur_coofile
                if markcoords and cooass: self.mark_coords(coofile=cooass)
                # self.__shiftmode = False

            elif key=='p':
                self.display(self.previous_image, ds9_target=ds9_target, wait_time=wait_time)
                cooass = self.get_cur_coofile
                if markcoords and cooass: self.mark_coords(coofile=cooass)

            elif key=='P':
                self._imdeque = deque(filenames2list(image, forcelist=True))
                self.display(self.current_image, ds9_target=ds9_target, wait_time=wait_time)
                cooass = self.get_cur_coofile
                if markcoords and cooass: self.mark_coords(coofile=cooass)
                # self.__shiftmode = False

            elif key == 'q':
                break

def imcorrection(image, instrument, bias='@biaslist', biaswin='full', biasbin=None, dark='@darklist', darkwin='full', darkbin=None, flat='@flatlist', flatwin='full', flatbin=None, illum='', fringe='',  output='__Auto__', flatnorm=True, configpath=None):
    '''
    perfomrs reduction: bias correction, flat fielding etc. For details go through the IRAF ccdproc document.
    :param image: image name(s): str, str with '*' or '@', list, tuple
    :param par: gain,readnoise: list or tuple
    :param bias: bias image(s): str, str with '*' or '@', list, tuple
    :param biaswin: If target image has a smaller window size than the bias image mention the window size
    :param biasbin: If target image is binned, bin: int
    :param dark : dark image(s): str, str with '*' or '@', list, tuple
    :param darkwin: same as biaswin
    :param darkbin: same as biasbin
    :param flat: flat image(s): str, str with '*' or '@', list, tuple
    :param flatwin: same as biaswin
    :param flatbin: same as biasbin
    :param output: new reduced image name(s): str, str with '*' or '@', list, tuple
    :param flatnorm: if set flat is normalized: bool
    :return: returns the output name(s), masterbias name, master dark name, master flat name, illum file name, fring file name
    '''
    parser = ConfigParser()
    if not configpath:
        configpath = os.path.join(_configpath, instrument)
    else:
        configpath = os.path.join(configpath, instrument)
    parser.read(configpath)
    try:
        gain = parser.get('values', 'gain')
    except:
        gain = parser.get('header_keys', 'gain')
    try:
        readnoise = parser.get('values', 'readnoise')
    except:
        readnoise = parser.get('header_keys', 'readnoise')
    par = [gain, readnoise]
    if output == '__Auto__':
        pref = ''
        if bias != '' and bias is not None:  pref += 'b'
        if dark != '' and dark is not None: pref += 'd'
        if flat != '' and flat is not None: pref += 'f'
        if illum != '' and illum is not None: pref += 'i'
        if fringe != '' and fringe is not None: pref += 'p'
        pref+='_'
        if type(image)==str and image[0]=='@': output = addpresuff2filename(image,pref)
        else: output = addpresuff2filename(filenames2list(image),pref)
        if len(filenames2list(bias, forcelist=True)) > 1: bias = master_bias_gen(par, bias)
        if len(filenames2list(bias, forcelist=True)) >1: dark = master_dark_gen(par, dark)
        if len(filenames2list(flat)) > 1: flat = master_flat_gen(par, flat, normalize=flatnorm)
        if biasbin is not None: winbin_fits(bias, win=biaswin, bin=biasbin, bintype=np.mean)
        if darkbin is not None: winbin_fits(dark, win=darkwin, bin=darkbin, bintype=np.mean)
        if flatbin is not None: winbin_fits(flat, win=flatwin, bin=flatbin, bintype=np.mean)
        Pyccdproc(image, bias, dark, flat, illum, fringe, output=output)
        return output, bias, dark, flat, illum, fringe

def alignimages(image, refimage=0, refcoord=0, coordfile='', coordsprefix='', ibox=200, fbox=5, shiftfile='', iterate=5, output='', verbose=True):
    '''
    Aligns the images, similar to IRAF imalign task, but easier and automatic
    :param image: image name(s): str, str with '*' or '@', list, tuple
    :param refimage: reference image or its index in the image list w.r.t. which the rest will be aligned: str ot int
    :param refcoord: reference co-ordinate index in the co-ordinates list w.rt. which the frames wil be aligned
    :param coordfile (optional): File containing the co-ordinates of the ref image. If skipped, the function will look for coordsprefix
    parameter: str
    :param coordsprefix: If coordfile is not skipped, prefix can be given here and the function will look for a file with name in the pattern
    prefix+refimname+'.coo', default value is '': str
    :param ibox: Size of box given, over which the farmes will be aligned: int
    :param fbox: Size of final box. It will check whether the misalignment after the process is less than this box size: int
    :param shiftfile (optional): File containing initial shifts of the images w.r.t. the ref image, if skipped, it will calculate the shifts: str
    :param iterate: Max iteration for alignment: int
    :param output: Output file name(s)
    :return: output file(s), coordinate file(s), maximum shift among the images
    '''
    inlist = filenames2list(image,forcelist=True)
    oplist = filenames2list(output,forcelist=True)
    if type(refimage)==int: refimage,irefimage = inlist[refimage], refimage
    else: refimage,irefimage = refimage, inlist.index(refimage)
    if coordfile=='': coordfile = get_coords_filename_from_prefix(image=refimage, prefix=coordsprefix)
    xref, yref = np.loadtxt(coordfile, usecols=(0, 1), unpack=False)[refcoord]
    if shiftfile=='':
        if os.path.exists('shiftlist'): os.remove('shiftlist')
        coo = get_psfcentroid(image,xref,yref,ibox)
        np.savetxt('shiftlist', coo-coo[irefimage])
        shiftfile = 'shiftlist'
    task = ir.imalign
    task.shifts = shiftfile
    task.shiftimages = 'yes'
    task.trimimages = 'no'
    task.interp_type = 'linear'
    task.niterate = iterate
    for op in filenames2list(output, forcelist=True):
        if os.path.exists(op) and op not in inlist: os.remove(op)
    if verbose:
        print('reference image:',refimage)
        print('reference co-ordinates:',(xref,yref))
        print('shiftfile:',shiftfile)
    task(input=image, reference=refimage, coords=coordfile, output=output, Stdout=os.devnull)
    coo = get_psfcentroid(output, xref, yref, fbox)
    maxshift = np.max(np.abs(coo-coo[irefimage]))
    return {'outputs': oplist, 'reference_coordinates_file': coordfile, 'maximum_misalignment': maxshift}

def save_matched_coords(image, coordfile='',output='',box=50,**kwargs):
    '''
    Finds and stores co-ordinates of the objects on all the frames based on the reference co-ordinates
    '''
    im = filenames2list(image)
    imlist = [im] if type(im)==str else im
    if output!='':
        output = filenames2list(output)
    urefcoo = True
    coordsprefix = kwargs.get('coordsprefix','')
    if coordfile=='':
        refim = kwargs.get('refimage', 0)
        if type(refim) == int:
            refim, irefim = imlist[refim], refim
        else:
            refim, irefim = refim, imlist.index(refim)
        coordfile = get_coords_filename_from_prefix(image=refim, prefix=coordsprefix)
        urefcoo = kwargs.get('update_refim_coords',True)
    if type(im)==str:
        if output=='':
            output = get_coords_filename_from_prefix(image=im,prefix=coordsprefix)
        if coordfile!=output and not urefcoo:
            os.system(f'cp {coordfile} {output}')
        if urefcoo:
            task = ir.center
            task.centerpars.calgorithm = 'centroid'
            task.centerpars.cbox = box
            task.centerpars.cmaxiter = 10
            task.centerpars.maxshift = 2
            task.plotfile = ''
            task.interactive = 'no'
            task.radplots = 'no'
            task.verify = 'no'
            task.update = 'no'
            task.verbose = 'no'
            if os.path.exists(j2s('__coo__')):
                os.remove(j2s('__coo__'))
            task(image=im, coords=coordfile, output=j2s('__coo__'))
            ir.txdump(textfiles=j2s('__coo__'), fields='xcenter,ycenter', expr='yes', headers='no', parameters='no', Stdout=output)
        return output
    if type(output)==str and output!='':
        warn('For multiple input images you need to give multiple output files either as a tuple or list of textfiles '
             'or a comma separated string of textfilenames or a list of filenames. Output file names will be set '
             'automatically for now.')
        output = ''
    oplist = []
    for i,im in enumerate(imlist):
        if output!='':
            j = list(range(len(imlist))).remove(irefim)[i]
            op = output[j]
        else:
            op = ''
        oplist.append(save_matched_coords(im,coordfile,op,box,**kwargs))
    return oplist







