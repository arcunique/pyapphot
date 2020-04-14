from setuptools import setup
import pyapphot
setup(name='pyapphot',
      version=pyapphot.__version__,
      description='Python package for photometric data reduction and aperture photometry',
      long_description=open('README.md').read(),
      classifiers=['Development Status :: 0 - Alpha',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3.7'],
      url='https://github.com/arcunique/pyapphot',
      author='Aritra Chakrabarty',
      author_email='arcunique1993@gmail.com',
      install_requires=['numpy', 'pyraf', 'astropy', 'matplotlib', 'imexam', 'pickle'],
      zip_safe=False)