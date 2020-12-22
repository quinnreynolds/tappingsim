from setuptools import setup
from pathlib import Path

here = Path(__file__).parent.resolve()

def read(rel_path):
    with open(here / rel_path, 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError('Unable to find version string.')

setup(name='tappingsim',
      version=get_version('tappingsim/__init__.py'),
      description="""Semi-empirical simulator for furnace tapping systems""",
      url='https://github.com/quinnreynolds/tappingsim',
      author='Quinn Reynolds',
      author_email='quinnr@mintek.co.za',
      license='Proprietary',
      packages=['tappingsim'],
      install_requires=['numpy', 'scipy'],
      zip_safe=False)
