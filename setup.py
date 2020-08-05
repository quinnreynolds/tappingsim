from setuptools import setup
from ladlesim import __version__

setup(name='ladlesim',
      version=__version__,
      description="""Semi-empirical simulator for furnace tapping systems""",
      url='https://github.com/quinnreynolds/ladlesim',
      author='Quinn Reynolds',
      author_email='quinnr@mintek.co.za',
      license='Proprietary',
      packages=['ladlesim'],
      install_requires=['numpy', 'scipy'],
      zip_safe=False)
