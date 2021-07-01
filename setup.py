from setuptools import setup, find_packages

setup(
    name='alminer',
    version='0.1.0',
    packages=find_packages(include=['alminer', 'alminer.*']),
    url='https://github.com/emerge-erc/ALminer',
    license='MIT',
    author='Aida Ahmadi',
    author_email='aahmadi@strw.leidenuniv.nl',
    description='ALminer: ALMA archive mining and visualization toolkit',
    install_requires=['numpy>=1.15', 'pandas>1.0', 'matplotlib>=3.3.0', 'pyvo>=1.1',
                      'astropy>=3.1.2', 'astroquery @ git+https://github.com/astropy/astroquery'],
    python_requires='>=3.6'
)
