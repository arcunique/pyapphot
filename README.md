pyapphot
======

[![Build Status](https://img.shields.io/badge/release-0.1-orange)](https://github.com/arcunique/pyapphot)
[![Python 2.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-371/)

pyapphot is a collection of classes and functions defined to perform reduction and analyses for aperture phototmetry.
This is primarily based on the packages [pyraf](https://astroconda.readthedocs.io/en/latest/) 
and [imexam](https://github.com/spacetelescope/imexam). The __Imexamine__ class of this package allows interactive
selection and deselection of objects from the frames, auto-detection of objects on multiple frames based on one reference frame,
auto-alignment of the frames etc. This also allows all the PyRAF (hence, IRAF) and imexam tasks apart from the methods defined
in this class. 

The __starPSF__ class of this package allows extraction of all the PSF information from all the frames which can be stored into a
 Python Pickle file. This class allows users to set their criteria to filter out the frames/objects with poor PSF and helps
 make a decision regarding the trade-off between the frames to be discarded and the objects to be discarded. 
 
 The __aperture_phot__ class allows users to perform aperture photometry and store the results into a Pickle file. Differential 
 photometry can also be performed using this class.
 

Author
------
* Aritra Chakrabarty (IIA, Bangalore)

Requirements
------------
* python>3.6
* numpy
* astropy
* pyraf
* imexam=0.9.1
* pickle
* matplotlib 
* photutils
* PyQt5

Instructions on installation and use
------------------------------------
Presently, the code is only available on [Github](https://github.com/arcunique/pyapphot). Either download the code or
use the following line on terminal to install using pip:\
pip install git+https://github.com/arcunique/pyapphot  #installs from the current master on this repo.

You can import the modules and classes by:

from pyapphot import image_management\
from pyapphot.apphot import starPSF\
from pyapphot.apphot import aperture_phot

Documentation of this package is underway.Some example Jupyter notebooks can be found in this 
package which demonstrate how to use thse classes and functions. This package has already been 
used to perform differential transit photometry, results from which can be found 
[here](https://doi.org/10.3847/1538-3881/ab24dd).






