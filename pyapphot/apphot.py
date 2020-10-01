''' This module calculates PSF property and perfomrs aperture photometry '''

import numpy as np
from pyraf import iraf as ir
from .utils import *
from astropy.io import fits
import pickle as pkl
from .image_management import gaussfit_psf
from . import configpath
from copy import deepcopy as copy

__SIZE_STATS__ = ['Radius','FWHM','GFWHM','MFWHM']
__STAT2PSFSTAT__ = {'radius':'Radius','fwhm':'FWHM','mfwhm':'MFWHM'}
__STAT_KEYS__ = list(map(str.lower,__SIZE_STATS__))+['mag','peak','ellip','pa','skysig']

def sky_sigma_gen(image, coordfile, dist=50, skybox=5):
    coords = np.loadtxt(coordfile, usecols=(0,1),unpack=False,ndmin=2)
    relpos = ([dist,dist],[dist,-dist],[-dist,-dist],[-dist,dist])
    image = filenames2list(image)
    if type(image)==str:
        skysig = [np.nan for _ in range(coords.shape[0])]
        if '[' in image and image.endswith(']'):
            imname, sci = image.split('[')[0], int(image.split('[')[1].rstrip(']'))
        else:
            imname, sci = image, 0
        data = fits.getdata(imname,sci)
        for c,coo in enumerate(coords):
            sigvarcoo = []
            for rp in relpos:
                try: sigvarcoo.append(np.std(data[int(coo[1]+rp[1]-skybox):int(coo[1]+rp[1]+skybox),int(coo[0]+rp[0]-skybox):int(coo[0]+rp[0]+skybox)]))
                except: sigvarcoo.append(np.nan)
            skysig[c] = np.nanmedian(sigvarcoo)
        return np.array(skysig)
    skysig = []
    for im in image:
        skysig.append(sky_sigma_gen(im,coordfile,dist,skybox))
    return np.array(skysig)

class starPSF(object):
    ''' Creates a database of the PSF properties and perfomrs filtering on the images based on the given criteria '''
    from iraf import obsutil
    psfstat = {}

    def __init__(self, image=None, coordfile='', calculate=False, statkey='all', loadfrom='', saveto=''):
        '''

        :param image: image name(s): str, str with '*' or '@', list, tuple
        :param coordfile: Co-ordinate file(s): str, str with '*' or '@', list, tuple
        :param calculate: If True it will calculate the properties defined in 'stakey' parameter: bool
        :param statkey: The properties to be calculated. Use the 'get_all_statkeys' to get the stat key names available,
        or otherwise, 'all' can be given to calculate all the properties: str or list
        :param loadfrom: If calculate=False, it will look for this parameter. It is the pickle filename in which the stats
        are already store: str
        :param saveto: Save the stats to: str
        '''
        if image:
            self.image = filenames2list(image, forcelist=True)
            self._image = image
        if coordfile:
            self.coordfile = filenames2list(coordfile, forcelist=True)
            self._coordfile = coordfile
        if calculate:
            self.psfstat = self.calculate(statkey=statkey)
            # if self.psfstat:
            #     self.psfstat['image'] = self.image
            #     self.psfstat['coordfile'] = self.coordfile
            #     self.psfstat['Nobject'] = [lendata(coo) for coo in self.coordfile]
        elif loadfrom and os.path.exists(loadfrom):
            self.psfstat = pkl.load(open(loadfrom,'rb'))
        if saveto and self.psfstat: pkl.dump(self.psfstat, open(saveto,'wb'))

    def calculate(self, image=None, coordfile='', statkey='all', saveto=''):
        if not statkey: return
        if not image: image = self._image
        if not image: return
        if not coordfile: coordfile = self._coordfile
        coordfile = filenames2list(coordfile, forcelist=True)
        if np.all([not lendata(coo) for coo in coordfile]): return
        task = ir.psfmeas
        task.coords = 'markall'
        task.wcs = 'logical'
        task.display = 'No'
        size = list(set([__STAT2PSFSTAT__[sk] if sk in __STAT2PSFSTAT__ else sk for sk in statkey]).intersection(
            __SIZE_STATS__)) if statkey != 'all' else __SIZE_STATS__
        if statkey == 'all': statkey = __STAT_KEYS__
        if not size: size = ['FWHM']
        if size or set(['mag', 'ellip', 'pa']).intersection(statkey):
            task.graphcur = 'graph'
            if not os.path.exists('graph'): print("'q'", file=open('graph', 'w'))
        out = '/tmp/_psfmeas_'
        image = filenames2list(image, forcelist=True)
        coordfile = filenames2list(coordfile, forcelist=True)
        if len(coordfile) == 1:
            coordfile *= len(image)
            Nrow = [lendata(coordfile[0])] * len(image)
        else:
            Nrow = [lendata(coordfile[i]) for i in range(len(image))]
        # if np.all(Nrow[i] == Nrow[0] for i in range(1, len(Nrow))):
        #     stat = {key: np.full((len(Nrow), Nrow[0]), np.nan) for key in __STAT_KEYS__}
        # else:
        stat = {key: np.array([np.full((Nrow[i],), np.nan) for i in range(len(Nrow))]) for key in __STAT_KEYS__}
        for i, im in enumerate(image):
            if not lendata(coordfile[i]): continue
            task.size = size[0]
            if size or set(['mag', 'ellip', 'pa']).intersection(statkey):
                task.imagecur = coordfile[i]
                skipim = False
                if os.path.exists(out): os.remove(out)
                task(images=im, Stdout=out)
                try:
                    with open(out) as fr:
                        lines = fr.readlines()[3:3 + Nrow[i]]
                    for key, val in zip(['mag', size[0].lower(), 'ellip', 'pa'], lines[0].split()[3:]):
                        stat[key][i, 0] = float(val)
                    if Nrow[i] > 1:
                        for l in range(1, Nrow[i]):
                            for key, val in zip(['mag', size[0].lower(), 'ellip', 'pa'], lines[l].split()[2:]):
                                stat[key][i, l] = float(val)
                except:
                    skipim = True
                    print('problem with file ' + im + ', storing Nan instead.')
                if len(size) > 1:
                    for s in range(1, len(size)):
                        if not skipim:
                            task.size = size[s]
                            if os.path.exists(out): os.remove(out)
                            task(images=im, Stdout=out)
                            try:
                                lines = open(out).readlines()[3:3 + Nrow[i]]
                                stat[size[s].lower()][i, 0] = float(lines[0].split()[4])
                                if Nrow[i] > 1:
                                    for l in range(1, Nrow[i]): stat[size[s].lower()][i, l] = float(lines[l].split()[3])
                            except:
                                skipim = True
                                print('problem with file ' + im + ', storing Nan instead.')
            if 'peak' in statkey:
                coo = np.loadtxt(coordfile[i], unpack=False)
                for c, (x, y) in enumerate(coo):
                    res = gaussfit_psf(im, x, y, )
                    if res: stat['peak'][i, c] = res[0]
            if 'skysig' in statkey:
                stat['skysig'][i] = sky_sigma_gen(im, coordfile[i], )
        if os.path.exists(out): os.remove(out)
        psfstat = {k: stat[k] for k in statkey if k in stat}
        psfstat['image'] = image
        psfstat['coordfile'] = coordfile
        psfstat['Nobject'] = [lendata(coo) for coo in coordfile]
        return psfstat

    def get_stat(self, image='all', obj='all', statkey='all', loadfrom='', return_dict=False, saveto=''):
        '''
        Get info about the selected stats for selected images from the loaded or already claculated stat databse
        '''
        if not statkey: return
        if not image: return
        if loadfrom and os.path.exists(loadfrom): _psfstat = pkl.load(open(loadfrom,'rb'))
        else: _psfstat = self.psfstat.copy()
        _image = _psfstat.pop('image')
        _coordfile = _psfstat.pop('coordfile')
        _Nobject = _psfstat.pop('Nobject')
        if image!='all': image = filenames2list(image, forcelist=True)
        else: image = _image
        imind = np.in1d(_image,image)
        imex = ','. join([im for im in image if im not in _image]) if sum(imind) else 'all input images'
        if imex:
            if loadfrom: print(f'No PSF stat found in {loadfrom} for {imex}')
            else: print(f'No PSF stats have been calculated for {imex}')
        if statkey=='all': statkey = _psfstat.keys()
        stat = {k: _psfstat[k][imind] for k in statkey if k in _psfstat}
        if obj!='all':
            if type(obj)==int: obj = [[obj] for _ in image]
            elif isinstance(obj[0], int): obj = [obj for _ in image]
            print(obj)
            # print(np.any(len(obj[i])!=len(obj[0]) for i in range(1, len(image))))
            for key in stat:
                # if stat[key].ndim==2 and np.any([len(obj[i])!=len(obj[0]) for i in range(1, len(image))]):
                stat[key] = [stat[key][i] for i in range(len(image))]
                for i in range(len(image)):
                    stat[key][i] = stat[key][i][list(obj[i])]
                stat[key] = np.array(stat[key])
        if return_dict or saveto:
            stat['image'] = image
            stat['coordfile'] = _coordfile
            stat['Nobject'] = [len(ob) for ob in obj]
        if saveto: pkl.dump(open(saveto,'wb'), stat)
        return stat if return_dict else [stat[key] for key in statkey]

    def select_imageNobj(self, criteria, objsynced=True, Nfiledelmax=None):
        '''
        Select images and objects which are going to be included for calculation depending on the criteria given
        :param criteria: a dictionary, keys are the stat keys and values are the lower and upper limits in tuple or list form: dict
        :param objsynced: Always needs to be True. Do not change it. Chnaging it will give no output but won't throw any error.
        :param Nfiledelmax: Maximum number of files that can be allowed to be deleted, accordingly the objects in the remaining files
        will be deleted, if the criteria are not met.
        :return: The images that met the criteria for each object and a dictionary of the images discarded and the objects discarded.
        '''
        psfstat = self.psfstat.copy()
        _image = psfstat.pop('image')
        _coordfile = psfstat.pop('coordfile')
        _Nobject = psfstat.pop('Nobject')
        for key in criteria:
            if key not in psfstat: criteria.pop(key)
        criteriamet = np.empty(len(criteria), dtype=object)
        for k,key in enumerate(criteria):
            if np.ndim(psfstat[key])==2: criteriamet[k] = (psfstat[key]>=criteria[key][0]) & (psfstat[key]<=criteria[key][1])
            else: criteriamet[k] = np.array([psfval>=criteria[key][0] & psfval<=criteria[key][1] for psfval in psfstat[key]])
        if objsynced:
            crinotmet = np.empty(_Nobject[0], dtype=object)
            immet = []
            if Nfiledelmax:
                imdel = np.empty(_Nobject[0], dtype=object)
                removeobj = []
            for i in range(_Nobject[0]):
                crinotmet[i] = np.any([np.logical_not(criteriamet[k][:,i]) for k in range(len(criteria))], axis=0)
                _image = np.array(_image)
                imnotmet = _image[crinotmet[i]]
                immet.append(list(_image[np.logical_not(crinotmet[i])]))
                if not np.any(crinotmet[i]): message = f'object {i}: conditions met by all the images.'
                elif np.all(crinotmet[i]): message = f'object {i}: conditions not met by any image.'
                else: message = f'object {i}: conditions not met by {len(imnotmet)} image(s).'
                print(message)
                if Nfiledelmax:
                    if len(imnotmet)<Nfiledelmax: imdel[i] = imnotmet
                    else:
                        imdel[i] = None
                        removeobj.append(i)
            if not Nfiledelmax: return immet, dict(delete_image_set=None, delete_image_array=None, remove_obj=None)
            imdelset = list(set([idl for idel in imdel if len(idel)!=0 for idl in idel]))
            if len(imdelset)<=Nfiledelmax: return immet, dict(delete_image_set=imdelset, delete_image_array=None, remove_obj=removeobj)
            return immet, dict(delete_image_array=immet, delete_image_set=None, remove_obj=removeobj)


def apphot(image,coordfile,instrument,fwhm,skysig,aperture,skypos,output,countrange=('',''),instrupath=None):
    from configparser import ConfigParser
    parser = ConfigParser()
    if not instrupath: parser.read(os.path.join(configpath,instrument))
    else: parser.read(os.path.join(instrupath,instrument))
    indefflag = {}
    try:
        gain = parser.get('header_keys', 'gain'), 'indef'
        gi = 0
    except:
        gain = '', parser.get('values', 'gain')
        gi = 1
    try: readnoise = parser.get('header_keys', 'readnoise'), 'indef'
    except: readnoise = '', parser.get('values', 'readnoise')
    try: exptime = parser.get('header_keys', 'exposure'), 'indef'
    except: exptime = '', parser.get('values', 'exposure')
    try: airmass = parser.get('header_keys', 'airmass'), 'indef'
    except: airmass = '', parser.get('values', 'airmass')
    try: obstime = parser.get('header_keys', 'obstime'), 'indef'
    except: obstime = '', parser.get('values', 'obstime')
    try: filter = parser.get('header_keys', 'filter'), 'indef'
    except: filter = '', parser.get('values', 'filter')
    skyapert = aperture[-1] + skypos[0], skypos[1]
    for var in ['gain','readnoise','exptime','airmass','obstime','filter']:
        if type(locals()[var][0])==str and locals()[var][0].lower() in ['','indef']\
                and type(locals()[var][1])==str and locals()[var][1] in ['','indef']: indefflag[var] = True
        else: indefflag[var] = False
    task = ir.phot
    task.coords = coordfile
    task.interactive = 'no'
    task.verify = 'no'
    task.update = 'no'
    task.verbose = 'no'
    task.cache = 'no'
    td = task.datapars
    td.scale = 1
    td.fwhmpsf = fwhm
    td.sigma = skysig
    td.emission = 'yes'
    td.datamin = countrange[0] if countrange[0] != '' else 'INDEF'
    td.datamax = countrange[1] if countrange[1] != '' else 'INDEF'
    td.noise = 'poisson'
    td.gain, td.epadu = gain
    td.ccdread, td.readnois, = readnoise
    td.exposure, td.itime = exptime
    td.airmass, td.xairmass = airmass
    td.filter, td.ifilter = filter
    td.obstime, td.otime = obstime
    task.fitskypars.annulus = skyapert[0]
    task.fitskypars.dannulus = skyapert[1]
    task.fitskypars.salgorithm = 'centroid'
    task.photpars.zmag = 25
    task.photpars.mkapert = 'no'
    task.photpars.apertures = ','.join(np.array(aperture, dtype=str))
    print(image)
    task(image=image, output=output)
    return indefflag, float(gain[gi])



class aperture_phot(object):
    ''' Performs aperture photometry, as well as differential photometry.'''
    _photdata = {}

    def __init__(self, image=None, coordfile=None, calculate=False, loadfrom='', saveto='', **kwargs):
        '''
        :param image: Image name(s): str, str with '*' or '@', list, tuple
        :param coordfile: Co-ordinate file(s):  str, str with '*' or '@', list, tuple
        :param calculate: if True it will perform photometry: bool
        :param loadfrom: If calculate==False, then this parameter is looked for, the pickle file name from which photometric
        results such as time, flux, flux-error, mag, mag-error : str
        :param saveto: pickle file to which photometric data to be stored : str
        '''
        if image is not None:
            self._image = copy(image)
            self.image = filenames2list(image, forcelist=True)
        if coordfile is not None:
            self._coordfile = coordfile
            self.coordfile = filenames2list(coordfile, forcelist=True)
        if calculate:
            self(self._image,self._coordfile,**kwargs)
        elif loadfrom:
            self.photdata = pkl.load(open(loadfrom,'rb'))
        if self.photdata and saveto:
            pkl.dump(self.photdata, open(saveto, 'wb'))

    def _getlatestfilename(self,filename):
        fname,fext = os.path.splitext(filename)
        files = filenames2list(fname+'*'+fext, forcelist=True)
        if files:
            num = int(files[-1].replace(fname,'').replace(fext,''))
            return fname+f'{num+1:03d}'+fext
        else:
            return fname+f'{1:03d}'+fext

    def __call__(self, image=[], coordfile='', saveto='', **kwargs):
        if 'imstat' not in kwargs: return
        _image = copy(image)
        image = filenames2list(image, forcelist=True)
        if not image:
            image = self.image
            _image = self._image
        _coordfile = copy(coordfile)
        coordfile = filenames2list(coordfile, forcelist=True)
        if not coordfile:
            coordfile = self.coordfile
            _coordfile = self._coordfile
        if not image or not coordfile: return
        output = kwargs.setdefault('output',self._getlatestfilename('photop.mag'))
        imstat = kwargs.pop('imstat')
        if type(imstat)==str: imstat = pkl.load(open(imstat,'rb'))
        imstat = filterdict(imstat, 'image', image)
        if type(imstat) == str: imstat = pkl.load(open(imstat, 'rb'))
        _ = kwargs.setdefault('fwhm',imstat['fwhm'].mean())
        _ = kwargs.setdefault('skysig',imstat['skysig'].mean())
        flag, gain = apphot(_image,_coordfile,**kwargs)
        output = filenames2list(output,forcelist=True)
        Nobj = imstat['Nobject']
        self._uniform = all([nobj==Nobj[0] for nobj in Nobj[1:]])
        if flag['obstime']: self.t = []
        if len(output)==1:
            res = list(self._getfrommagfile(output[0], get_t=not flag['obstime'], get_apert=False, get_others=True))
            if not flag['obstime']:
                t = res.pop(0)
                self.t = np.array([t[i] for i in np.insert(np.cumsum(Nobj[:-1]),0,0)])
            exptime, stdev, nsky, flux, area, mag = res
            exptime = segregate(exptime, Nobj)
            stdev = segregate(stdev, Nobj)
            nsky = segregate(nsky, Nobj)
            self.flux = segregate(flux, Nobj)
            area = segregate(area, Nobj)
            self.mag = segregate(mag, Nobj)
        else:
            exptime,stdev,nsky,area = [[] for _ in range(4)]
            self.flux, self.mag = [], []
            if not flag['obstime']: self.t = np.array([])
            for op in output:
                res = list(self._getfrommagfile(op, get_t=not flag['obstime'], get_apert=False, get_others=True))
                if not flag['obstime']: self.t = np.append(self.t, res.pop(0)[0])
                exptime.append(res[0])
                stdev.append(res[1])
                nsky.append(res[2])
                self.flux.append(res[3])
                area.append(res[4])
                self.mag.append(res[5])
        self.ferr = [None for _ in self.flux]
        self.merr = [None for _ in self.flux]
        for i in range(len(self.flux)):
            self.flux[i] = self.flux[i] / exptime[i][:,None]
            self.ferr[i] = np.transpose(np.sqrt(self.flux[i].T/gain + area[i].T * stdev[i]**2 + area[i].T**2 * stdev[i]**2 / nsky[i]))
            self.merr [i] = 1.0857 * self.ferr[i] / flux[i]
        if self._uniform:
            self.flux = np.array(self.flux)
            self.ferr = np.array(self.ferr)
            self.mag = np.array(self.mag)
            self.merr = np.array(self.merr)
        apert = self._getfrommagfile(output[0],get_t=False,get_apert=True,get_others=False)[0]
        self._photdata = dict(t=self.t,flux=self.flux,ferr=self.ferr,mag=self.mag,merr=self.merr,image=image,Nobject=Nobj,aperture=apert,diminfo={'t':0,'image':0,'object':1,'aperture':2})
        if saveto: pkl.dump(self._photdata, open(saveto, 'wb'))

    @staticmethod
    def _getfrommagfile(magfile,get_t=True,get_apert=True,get_others=True):
        res = []
        if get_t:
            ir.txdump(magfile, fields='otime', expr='yes', headers=False, parameters=False, Stdout=j2s('__obs__'))
            try: t = np.loadtxt(j2s('__obs__'),)
            except ValueError: t = np.loadtxt(j2s('__obs__'), usecols=0, dtype=object)
            res.append(t)
            os.remove(j2s('__obs__'))
        if get_others:
            print(magfile)
            ir.txdump(magfile, fields='itime,stdev,nsky', expr='yes', headers=False, parameters=False, Stdout=j2s('__isn__'))
            exptime, stdev, nsky = np.loadtxt(j2s('__isn__'), unpack=True)
            os.remove(j2s('__isn__'))
            ir.txdump(magfile, fields='flux', expr='yes', headers=False, parameters=False, Stdout=j2s('__flux__'))
            flux = np.loadtxt(j2s('__flux__'), ndmin=2)
            os.remove(j2s('__flux__'))
            ir.txdump(magfile, fields='area', expr='yes', headers=False, parameters=False, Stdout=j2s('__area__'))
            area = np.loadtxt(j2s('__area__'), ndmin=2)
            os.remove(j2s('__area__'))
            ir.txdump(magfile, fields='mag', expr='yes', headers=False, parameters=False, Stdout=j2s('__mag__'))
            mag = np.loadtxt(j2s('__mag__'), ndmin=2)
            os.remove(j2s('__mag__'))
            res += [exptime, stdev, nsky,flux,area,mag]
        if get_apert:
            ir.txdump(magfile, fields='rapert', expr='yes', headers=False, parameters=False, Stdout=j2s('__apert__'))
            apert = np.loadtxt(j2s('__apert__'))[0]
            res.append(apert)
            os.remove(j2s('__apert__'))
        return res

    @property
    def photdata(self): return self._photdata

    @photdata.setter
    def photdata(self, data):
        self._photdata.update(data)
        for key in data:
            self.__setattr__(key, data.get(key, []))

    def differential_photometry(self,image=None,objpair=(0,1),aperture=5):
        '''
        Perfomrs differential photometry
        :param image: Image(s): str, str with '*' or '@', list, tuple
        :param objpair: Object index pair for which differential photometry is to be performed. It will claculate
        flux[objpair[0]]/flux[objpair[1]] and mag[objpar[0]-mag[objpair[1]].
        :param aperture: The aperture for which the calculation to be done
        :return: time, differential flux, fux-error, differential mag, mag-error
        '''
        image = filenames2list(image)
        if not image:
            photdata = self._photdata.copy()
        else:
            photdata = filterdict(self.photdata,'image',image)
        photdata = filterdict(photdata,'aperture',aperture)
        flux, ferr = photdata['flux'], photdata['ferr']
        mag, merr = photdata['mag'], photdata['merr']
        dflux = flux[:,objpair[0],:]/flux[:,objpair[1],:]
        dferr = dflux*np.sqrt((ferr[:,objpair[0],:]/flux[:,objpair[0],:])**2+(ferr[:,objpair[1],:]/flux[:,objpair[1],:])**2)
        dmag = mag[:,objpair[0],:]-mag[:,objpair[1],:]
        dmerr = np.sqrt(merr[:,objpair[0],:]**2+merr[:,objpair[1],:]**2)
        t = photdata['t']
        return t,dflux,dferr,dmag,dmerr







