from setuptools import setup
import os

long_description = """
A lightweight Python version of Liger for calculating allowed dwell times given spacecraft thermal constraints.

* Uses a Monte-Carlo approach to sample a realistic ensemble of perigee-exit starting temperatures.
* Generates both typically achievable (50%) and best-case (90%) dwell times.
* Includes PSMC, Minus-Z and PLINE models.
* Constraints implemented as Python classes derived from a base ConstraintModel class.
"""

__version__ = open(os.path.join(os.path.dirname(__file__), 'pyger/VERSION')).read().strip()

setup(name='pyger',
      version=__version__,
      description='Calculate Chandra dwell times given thermal constraints',
      long_description=long_description,
      author='Tom Aldcroft',
      author_email='aldcroft@head.cfa.harvard.edu',
      url='http://cxc.harvard.edu/mta/ASPECT/tool_doc/pyger',
      download_url='',
      license='BSD',
      platforms=['any'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Programming Language :: Python :: 2',
          ],
      packages=['pyger'],
      include_package_data=True,
      package_data = {'': ['*.dat', '*.json', 'VERSION']},
      )
