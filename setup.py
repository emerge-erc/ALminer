from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='alminer',
    version='0.1.0',
    author='Aida Ahmadi',
    author_email='aahmadi@strw.leidenuniv.nl',
    description='ALminer: ALMA archive mining and visualization toolkit',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(include=['alminer', 'alminer.*']),
    url='https://github.com/emerge-erc/ALminer',
    project_urls={
        "Bug Tracker": "https://github.com/emerge-erc/ALminer/issues"
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license='MIT',
    install_requires=['numpy>=1.15', 'pandas>1.0', 'matplotlib>=3.3.0', 'pyvo>=1.1',
                      'astropy>=3.1.2', 'astroquery @ git+https://github.com/astropy/astroquery'],
    python_requires='>=3.6'
)
