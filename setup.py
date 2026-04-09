#!/usr/bin/env python3
"""
PhyloSOLID - Tree building from single-cell sequencing data
"""

from setuptools import setup, find_packages
import os

# Read version
version_file = os.path.join(os.path.dirname(__file__), 'src', '__version__.py')
with open(version_file, 'r') as f:
    exec(f.read())

# Read README
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open('requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='phylosolid',
    version=__version__,
    description='Tree building from single-cell sequencing data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/douymLab/PhyloSOLID',
    
    packages=find_packages(),  # 注意：这里改回 find_packages() 包含 cli 和 src
    package_dir={'': '.'},      # 项目根目录
    
    include_package_data=True,
    package_data={
        '': ['scripts/**/*.sh', 'scripts/**/*.R', 'scripts/**/*.py', 'config/*.yaml', 'config/*.template'],
        'src': ['**/*.py'],     # 确保包含 run_phylosilid_fullTree_binput.py
    },
    
    python_requires='>=3.7',
    install_requires=requirements,
    
    entry_points={
        'console_scripts': [
            'phylosolid=cli.main:main',
        ],
    },
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    
    license='MIT',
)