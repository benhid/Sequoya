from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='pyM2SA',
    version='0.3.0',
    description='Solving Multiple Sequence Alignments with Python',
    url='https://github.com/benhid/pyM2SA',
    author='Antonio Benítez-Hidalgo',
    author_email='antonio.b@uma.es',
    license='MIT',
    python_requires='>=3',
    packages=find_packages(exclude=["test.*", "tests"]),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6'
    ],
    install_requires=[
        'jmetalpy==0.5.0',
        'pymsa==1.0.0'
    ]
)