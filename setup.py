from setuptools import setup, find_packages

setup(
    name='summer_py',
    version='1.1.0',
    packages=['summer_py',],
	url='https://github.com/jtrauer/summer_py',
    license='MIT',
    author='James Trauer',
    author_email='james.trauer@monash.edu',
    install_requires=['scipy==1.1.0',
                      'graphviz>=0.4.10',
                      'SQLAlchemy>=1.1.18',],
	description='General structure for creating epidemiological models in R'
)
