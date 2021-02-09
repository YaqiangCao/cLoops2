from glob import glob
from setuptools import setup, find_packages
#from Cython.Build import cythonize
from cLoops2.cLoops2 import __version__

scripts = glob("scripts/*.py")

setup(
    name='cLoops2',
    version=__version__,
    author='Yaqiang Cao',
    author_email='caoyaqiang0410@gmail.com',
    url='https://github.com/YaqiangCao/cLoops2',
    keywords='peak-calling loop-calling Hi-Trac interaction visualization',
    description='Loop-calling and peak-calling for sequencing-based interaction data, including related analysis utilities.',
    classifiers=[
        'Environment :: Console',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    project_urls={
                    'Source': 'https://github.com/YaqiangCao/cLoops2',
    },
    packages=find_packages(exclude=['tests', 'docs',"example"]),
    long_description=open('cLoops2_pip_readme.md').read(),
    long_description_content_type="text/markdown",
    setup_requires=["joblib", "numpy", "seaborn", "pandas", "scipy","scikit-learn","matplotlib","tqdm","pyBigWig","networkx"],
    install_requires=["joblib", "numpy", "seaborn", "pandas", "scipy","scikit-learn","matplotlib","tqdm","pyBigWig","networkx"],
    entry_points={
        'console_scripts': [
            'cLoops2=cLoops2.cLoops2:main',
        ],
    },
    scripts = scripts,
    python_requires='>=3',
)
