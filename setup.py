from glob import glob
from setuptools import setup, find_packages
#from Cython.Build import cythonize
from cLoops2.cLoops2 import __version__

scripts = glob("scripts/*.py")

setup(
    name='cLoops2',
    version=__version__,
    author=['Yaqiang Cao'],
    author_email=['caoyaqiang0410@gmail.com'],
    url='https://github.com/YaqiangCao/cLoops2',
    description='Loop-calling and peak-calling for sequencing-based interaction data, including related analysis utilities',
    classifiers=[
        'Environment :: Console',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages(exclude=['tests', 'docs']),
    long_description=open('README.md').read(),
    setup_requires=["joblib", "numpy", "seaborn", "pandas", "scipy","scikit-learn","matplotlib","tqdm","pyBigWig"],
    entry_points={
        'console_scripts': [
            'cLoops2=cLoops2.cLoops2:main',
        ],
    },
    scripts = scripts,
)
