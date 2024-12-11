from setuptools import setup, find_packages

setup(
    name='ScanyCraft',
    version='0.1',
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
