from setuptools import setup, find_packages

setup(
    name='ScanyCraft',
    version='0.1',
    description="A package for [your description here]",
    author="Hanzhou Wang",
    author_email="wangh5@uthscsa.edu",
    url="https://github.com/whzemuch/ScanyCraft",
    packages=find_packages(),
    install_requires=[
        'numpy',        # Specify minimum version
        'scipy',         # Specify minimum version
        'pandas',        # Specify minimum version
        'matplotlib>=3.4.0',    # Specify minimum version
        'scanpy>=1.8.0',        # Specify minimum version
        'scikit-posthocs', # Specify minimum version
        'tqdm',                  # No version specified, will install the latest
        'matplotlib-venn>=0.11.6' # Specify minimum version
    ],
)
