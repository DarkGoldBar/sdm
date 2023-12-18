from setuptools import setup, find_packages
import sdm

setup(
    name='sdm',
    version=sdm.__version__,
    author='bochen.li',
    author_email='bochen.li@xtalpi.com',
    description='San-D Molecule (SDM)',
    license='MIT',
    scripts=['scripts/sdm'],
    packages=find_packages(),
    install_requires=['spglib', 'numpy', 'networkx', 'pyyaml', 'jinja2'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 2 - Developing',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry']
)
