from setuptools import setup, find_packages

setup(
    name='rosesim',
    version='0.1.1',
    description='Image simulations for semi-resolved dwarf galaxies for Roman',
    author='Jiaxuan Li',
    author_email='jiaxuanl@princeton.edu',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy',
        'astroquery',
        'artpop @ git+https://github.com/AstroJacobLi/ArtPop.git',
        'asdf',
        'astrocut',
        'bs4',
        'requests',
        'tqdm',
        'fire',
        'sep',
        'roman_datamodels',
        'romanisim @ git+https://github.com/AstroJacobLi/romanisim.git',
    ],
    entry_points={
        'console_scripts': [
            'rosesim_sky=rosesim.scripts.sim_sky:main',
            'rosesim_gal=rosesim.scripts.sim_gal:main',
        ]
    },
    python_requires='>=3.7',
    tests_require=['pytest'],
)
