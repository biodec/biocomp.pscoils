from setuptools import setup, find_packages
version = '0.1'

setup(name='biocomp.pscoils',
      version=version,
      description="PS-COILS",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Bio',
      author='Piero Fariselli',
      author_email='piero@biocomp.unibo.it',
      url='http://www.biocomp.unibo.it/piero/PS-COILS/',
      license='GPLv3',
      packages=find_packages('src',exclude=['ez_setup', 'examples', 'tests']),
      package_dir={'':'src'},
      namespace_packages=['biocomp'],
      scripts=['scripts/pscoils',],
      test_suite='tests',
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
