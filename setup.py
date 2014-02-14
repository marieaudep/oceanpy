try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name = 'oceanpy',
      version = '0.1.0',
      author = 'Jonas Bluethgen',
      author_email = 'bluthgen@nbi.ku.dk',
      packages = ['oceanpy'],
      url = 'http://www.gfy.ku.dk/~bluthgen',
      license = 'LICENSE.txt',
      description = 'Physical Oceanography for Python',
      long_description = open('README.md').read(),
      )
