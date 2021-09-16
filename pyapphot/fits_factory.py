from astropy.io import fits
from astropy.time import Time
import os
from configparser import ConfigParser
from datetime import datetime, timedelta as delta
import numpy as np
from .utils import filenames2list
from . import _configpath

def fits_details(files='*.fits',keys='all', show=True, saveto=None):
    files = filenames2list(files, forcelist=True)
    if show:
        print('# ' + ' '.join(keys) if keys!='all' else '# all headers')
    if saveto is not None:
        infolist = ['# ' + ' '.join(keys)] if keys!='all' else '# all headers'
    for i, fn in enumerate(files):
        info = []
        if keys is 'all':
            if i==0:
                keys = fits.getheader(fn).keys()
        for key in keys:
            if key.lower()=='name':
                info.append(os.path.basename(fn))
            else:
                info.append(str(fits.getheader(fn)[key]))
        if show:
            print(' '.join(info))
    if saveto is not None:
        infolist.append(' '.join(info))
        with open(saveto,'w') as f:
            for inf in infolist:
                print(inf, file=f)


def update_fits_header(files, modify={}, remove=()):
    files = filenames2list(files, forcelist=True)
    for file in files:
        data, head = fits.getdata(file, header=True)
        for key in modify:
            if not isinstance(modify[key], dict):
                head[key] = modify[key]
            else:
                func = modify[key].get('func', None)
                args = modify[key].get('args', None)
                if func and args is None:
                    head[key] = func()
                if func is None and isinstance(args, str):
                    head[key] = head[args]
                if func is not None and args is not None:
                    if type(args) == str:
                        head[key] = func(head[args])
                    else:
                        head[key] = func(*[head[arg] for arg in args])
        if remove:
            for key in remove:
                if key in head:
                    head.remove(key)
        fits.update(file, data, header=head)

def _updateDatefromUT(date,ut):
    t = Time(date)
    ut = ':'.join(ut.split(':')[:2] + [str(np.round(float(ut.split(':')[2]), 4))])
    return str((t.datetime + delta(hours=-(datetime.strptime(ut, '%H:%M:%S.%f').hour-t.datetime.hour))).date())+'T'+ut

def _change_date_sep(date, sep):
    return '-'.join(date.split(sep))


class instrument_specific_fits_update:

    def __init__(self, instrument, configpath=None):
        self.instru = instrument
        self.parser = ConfigParser()
        if not configpath:
            configpath = os.path.join(_configpath, instrument)
        else:
            configpath = os.path.join(configpath, instrument)
        self.parser.read(configpath)

    def update_dateobs_with_correct_ut(self, files='*.fits'):
        files = filenames2list(files, forcelist=True)
        update_fits_header(files, modify={'date-obs': {'func': _updateDatefromUT, 'args': ('date-obs', 'ut')}})

    def setjd(self, files='*.fits'):
        func = lambda date: str(Time(date, format='isot', scale='utc').jd)
        update_fits_header(files, modify={'date-obs': {'func': func, 'args': 'date-obs'}})

    def setisot(self, files='*.fits'):
        func = lambda date: str(Time(date, format='jd', scale='utc').isot)
        update_fits_header(files, modify={'date-obs': {'func': func, 'args': 'date-obs'}})

    def concat_date_ut(self, files='*.fits'):
        func = lambda date, ut: date+'T'+str(delta(seconds=ut))
        update_fits_header(files, modify={'date-obs': {'func': func, 'args': ('date-obs', 'ut')}})

    def add_extension(self, files, extn='fits'):
        files = filenames2list(files, forcelist=True)
        for file in files:
            os.rename(file, file + '.'+extn)

    def segregate_imagetypes(self, files, objectnames=None, objectext=None):
        files = filenames2list(files, forcelist=True)
        biaslist = []
        flatdict = {}
        objectdict = {}
        objkw = self.parser.get('header_keys', 'object')
        filtkw = self.parser.get('header_keys', 'filter')
        for file in files:
            h = fits.getheader(file)
            obj = h[objkw]
            filt = h[filtkw]
            obj = obj.replace(' ', '_')
            filt = filt.replace(' ', '_')
            # print(file)
            # print('obj', obj)
            if 'bias' in obj.lower():
                biaslist.append(file)
            elif 'flat' in obj.lower():
                if filt not in flatdict:
                    flatdict[filt] = []
                flatdict[filt].append(file)
            elif objectnames is None or obj.lower() in objectnames:
                if obj not in objectdict:
                    objectdict[obj] = {}
                if filt not in objectdict[obj]:
                    objectdict[obj][filt] = []
                objectdict[obj][filt].append(file)
            # print(biaslist)
            # print(flatdict)
            # print(objectdict)
            if biaslist:
                with open('biaslist', 'w') as fw:
                    fw.writelines([file+'\n' for file in biaslist])
            if flatdict:
                for filt in flatdict:
                    with open('flatlist_'+filt, 'w') as fw:
                        fw.writelines([file + '\n' for file in flatdict[filt]])
            if objectdict:
                for obj in objectdict:
                    for filt in objectdict[obj]:
                        with open(obj+'_'+filt, 'w') as fw:
                            if type(objectext) == int:
                                file = file+'['+str(objectext)+']'
                            fw.writelines([file + '\n' for file in objectdict[obj][filt]])

    def covert_3dto2d(self, files):
        files = filenames2list(files, forcelist=True)
        for file in files:
            d, h = fits.getdata(file, header=True)
            d = d[0,:,:]
            h.remove('naxis3')
            fits.update(file, d, header=h)

    def execute_all(self, files, obsdate2jd=True):
        if self.instru=='JCBT_UKATC':
            self.update_dateobs_with_correct_ut(files)
            self.covert_3dto2d(files)
        if self.instru=='HFOSC':
            self.add_extension(files)
            update_fits_header(files, modify={'observat': 'iao'}, remove=('naxis3',))
            self.concat_date_ut(files)
        if self.instru=='HFOSC2':
            self.add_extension(files)
            update_fits_header(files, modify={'observat': 'iao'}, remove=('naxis3',))
        if self.instru=='TIRSPEC':
            update_fits_header(files, modify={'observat': 'iao'})
            _change_date_sep(files, '.')
        if obsdate2jd:
            self.setjd(files)


