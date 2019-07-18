from setuptools import setup, find_packages

setup(
    name='SUMMER',
    version='0.2',
    packages=find_packages(),
    url='https://github.com/jtrauer/summer',
    license='MIT',
    author='James Trauer',
    author_email='james.trauer@monash.edu',
    install_requires=['scipy>=1.1.0',
                      'graphviz>=0.4.10',
                      'SQLAlchemy>=1.1.18',
	              'pymc==2.3.7',	
                      'pandas==0.24.2',],
    description='General structure for creating epidemiological models in R'
)