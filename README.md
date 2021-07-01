<img src="images/ALminer_logo_header.jpg" alt="ALminer" align="center"/>

# ALminer: ALMA Archive Mining & Visualization Toolkit

`alminer` is a Python-based code to effectively query, analyse, and visualize the [ALMA science archive](https://almascience.eso.org/aq/). It also allows users to directly download ALMA data products and/or raw data for further image processing.

`alminer` has been developed through a collaboration between [Allegro](https://www.alma-allegro.nl/), the ALMA Regional Centre in The Netherlands, and the University of Vienna as part of the [EMERGE-StG project](https://emerge.alvarohacar.com). 


## Installation

The easiest way to install `alminer` is with `pip`:

```pip install alminer```

#### Dependencies

The dependencies are `numpy`, `pandas`, `matplotlib`, `pyvo`, [`astropy`](https://www.astropy.org/) version 3.1.2 or higher, and [`astroquery`](https://astroquery.readthedocs.io/en/latest/) version 0.4.2.dev6649 or higher. The `astroquery` package is only used for downloading data from the ALMA archive. The strict requirement to have its newest version is due to recent changes made to the ALMA archive. `alminer` works in Python 3.


## Usage

We have created an extensive [tutorial Jupyter Notebook](https://github.com/emerge-erc/ALminer/blob/main/notebooks/tutorial/ALminer_tutorial.ipynb) where all `alminer` features have been highlighted. This is an excellent starting points to get familiar with all the possibilities. We aim to keep adding new notebooks relevant for various sub-fields in the future. A glossery of all functions is provided at the bottom of this tutorial Jupyter Notebook.


## Acknowledgements

If you use `alminer` as part of your research, please consider adding the following ASCL citation to your work:

```
@article{TBD}

```

 `alminer` makes use of different routines in [Astropy](https://www.astropy.org/) and [Astroquery](https://astroquery.readthedocs.io/en/latest/). Please also consider citing the following papers:
- Astropy: [Astropy Collaboration et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A/abstract) <br>
- Astroquery: [Ginsburg et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....157...98G/abstract)

We also acknowledge the work of Leiden University M.Sc. students, Robin Mentel and David van Dop, who contributed to early versions of this work.

## Contact us

If you encounter issues, please [open an issue](https://github.com/emerge-erc/ALminer/issues). 

If you have suggestions for improvement or would like to collaborate with us on this project, please e-mail [Aida Ahmadi](mailto:aahmadi@strw.leidenuniv.nl) and [Alvaro Hacar](mailto:alvaro.hacar@univie.ac.at).

<img src="images/UniVie_logo.jpg" alt="University of Vienna" width= "350px"/><img src="images/ERC_logo.jpg" alt="ERC" width= "350px"/><img src="images/Allegro_logo.png" alt="Allegro"  width= "350px"/>
