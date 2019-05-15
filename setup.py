from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='Sequoya',
    version='0.9.0',
    description='Solving Multiple Sequence Alignments with Python',
    author='Antonio Benítez-Hidalgo',
    author_email='antonio.b@uma.es',
    maintainer='Antonio Benítez-Hidalgo',
    maintainer_email='antonio.b@uma.es',
    python_requires='>=3',
    license='MIT',
    url='https://github.com/benhid/Sequoya',
    long_description=open('README.md').read(),
    packages=find_packages(exclude=["test.*", "tests"]),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6'
    ],
    install_requires=[
        'jmetalpy==0.9.0',
        'pymsa==0.5.1',
        'dask[complete]==1.2.2',
        'bokeh==1.1.0'
    ]
)
