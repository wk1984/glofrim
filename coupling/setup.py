from setuptools import setup, find_packages
import os

setup(name='coupling_PCR_FM_2way',
      version='2.0',
      description="functions used for two-directional coupling between hydrology and hydrodynamics",
      long_desciption="""\
""",
      classifiers=[], 
      keywords='hydro',
      author='J.M. Hoch',
      author_email='jannis.hoch@deltares.nl',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ]
      )


