from setuptools import setup, find_packages

__version__ = '1.0.0'

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='switch',
    version=__version__,
    description='Hybrid python/R/bash package for IDG data engineering & data operations',
    url='https://github.com/teslajoy/gather-app',
    author='Nasim Sanati',
    author_email='nasim@plenary.org',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': ['switch = switch.__main__:main']
    },
    install_requires=[
        'requests',
    ],
    extras_require={
        'pandas': ["pandas==0.24.2"],
        'json': ["json5==0.8.4"],
        'reactome2py': ["reactome2py=0.0.3"],
        'scipy': ["scipy=1.3.0"],
        'statsmodels': ["statsmodels=0.10.0"],
    },
    tests_require=['pytest'],
    long_description=long_description,
    long_description_content_type='text/markdown',
)
