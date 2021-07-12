<img src="images/ALminer_logo_header.jpg" alt="ALminer" align="center"/>

# ALminer: ALMA Archive Mining & Visualization Toolkit

`alminer` is a Python-based code to effectively query, analyse, and visualize the [ALMA science archive](https://almascience.eso.org/aq/). It also allows users to directly download ALMA data products and/or raw data for further image processing.

## Installation

The easiest way to install `alminer` is with `pip`:

```pip install alminer```

To obtain the most recent version of the code from Github:

```pip install https://github.com/emerge-erc/ALminer/archive/refs/heads/main.zip```

Or clone and install from source:
```
# If you have a Github account:
git clone git@github.com:emerge-erc/ALminer.git
# If you do not:
git clone https://github.com/emerge-erc/ALminer.git

# After cloning:
cd ALminer
pip install .
```

Note that depending on your setup, you may need to use pip3.

### Dependencies

The dependencies are `numpy`, `matplotlib`, [`pandas`](https://pandas.pydata.org/), [`pyvo`](https://pyvo.readthedocs.io/en/latest/), [`astropy`](https://www.astropy.org/) version 3.1.2 or higher, and [`astroquery`](https://astroquery.readthedocs.io/en/latest/) version 0.4.2.dev6649 or higher. We only use the `astroquery` package for downloading data from the ALMA archive. The strict requirement to have its newest version is due to recent changes made to the ALMA archive. `alminer` works in Python 3. 


## Getting started

We have created an extensive [tutorial Jupyter Notebook](https://nbviewer.jupyter.org/github/emerge-erc/ALminer/blob/main/notebooks/tutorial/ALminer_tutorial.ipynb?flush_cache=True) where all `alminer` features have been highlighted. This is an excellent starting point to get familiar with all the possibilities; a glossery of all functions is provided at the bottom of this notebook. We highly recommend working in a [Jupyter notebook environment](https://jupyter.org/install) in order to make use of `alminer`'s visualization tools. We aim to keep adding new notebooks relevant for various sub-fields in the future. 

To work with the tutorial notebook interactively:

[![badge](https://img.shields.io/badge/launch-Jupyter%20Notebook-579ACA.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/emerge-erc/ALminer/main?filepath=notebooks/tutorial/ALminer_tutorial.ipynb)

## Acknowledgements

`alminer` has been developed through a collaboration between [Allegro](https://www.alma-allegro.nl/), the ALMA Regional Centre in The Netherlands, and the University of Vienna as part of the [EMERGE-StG project](https://emerge.alvarohacar.com). This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 851435).

If you use `alminer` as part of your research, please consider citing this [ASCL article](https://ascl.net/code/v/2971) (ADS reference will be added to the Github page when available).

 `alminer` makes use of different routines in [Astropy](https://www.astropy.org/) and [Astroquery](https://astroquery.readthedocs.io/en/latest/). Please also consider citing the following papers:
- Astropy: [Astropy Collaboration et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A/abstract) <br>
- Astroquery: [Ginsburg et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....157...98G/abstract)

We also acknowledge the work of Leiden University M.Sc. students, Robin Mentel and David van Dop, who contributed to early versions of this work. 

## Contact us

If you encounter issues, please [open an issue](https://github.com/emerge-erc/ALminer/issues). 

If you have suggestions for improvement or would like to collaborate with us on this project, please e-mail [Aida Ahmadi](mailto:aahmadi@strw.leidenuniv.nl) and [Alvaro Hacar](mailto:alvaro.hacar@univie.ac.at).

<img src="images/UniVie_logo.jpg" alt="University of Vienna" width= "280px" hspace="10px"/><img src="images/ERC_logo.jpg" alt="ERC" width= "200px" hspace="10px"/><img src="images/Allegro_logo.png" alt="Allegro"  width= "280px" hspace="10px"/>
