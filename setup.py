from setuptools import setup

setup(name='ladlesim',
      version='0.3',
      description="""Simple mass-balance simulator for tapping ladle filling""",
      url='https://kittychunk@bitbucket.org/kittychunk/ladlesim.git',
      author='Quinn Reynolds',
      author_email='quinnr@mintek.co.za',
      license='Proprietary',
      packages=['ladlesim'],
      install_requires=['numpy', 'scipy'],
      zip_safe=False)
