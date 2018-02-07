from setuptools import setup, find_packages

setup(
    name='pyM2SA',
    version='0.1.0',
    description='Solving Multiple Sequence Alignments with Python',
    url='https://github.com/benhid/pyM2SA',
    author='Antonio Benitez, Antonio J. Nebro',
    author_email='antonio.benitez@lcc.uma.es',
    license='MIT',
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6'
    ],
    packages=find_packages(exclude=["test.*", "tests"]),
)