.. image:: ../images/ALminer_logo_header.jpg

ALminer: ALMA Archive Mining & Visualization Toolkit
====================================================

``alminer`` is a Python-based code to effectively query, analyse, and
visualize the `ALMA science archive`_. It also allows users to
directly download ALMA data products and/or raw data for further image
processing.

Installation
------------

The easiest way to install ``alminer`` is with ``pip``:

``pip install alminer``

To obtain the most recent version of the code from `GitHub`_:

``pip install https://github.com/emerge-erc/ALminer/archive/refs/heads/main.zip``

Or clone and install from source:

::

    # If you have a Github account:
    git clone git@github.com:emerge-erc/ALminer.git
    # If you do not:
    git clone https://github.com/emerge-erc/ALminer.git

    # After cloning:
    cd ALminer
    pip install .

Note that depending on your setup, you may need to use pip3.

Dependencies
~~~~~~~~~~~~

The dependencies are ``numpy``, ``matplotlib``,
`pandas`_, `pyvo`_, `astropy`_ version 3.1.2 or higher, and
`astroquery`_ version 0.4.2.dev6649 or higher. We only use the ``astroquery`` package
for downloading data from the ALMA archive. The strict requirement to
have its newest version is due to recent changes made to the ALMA
archive. ``alminer`` works in Python 3.


Getting started
---------------

We have created an extensive `tutorial Jupyter Notebook`_
where all ``alminer`` features have been highlighted. This is an
excellent starting point to get familiar with all the possibilities; a
glossary of all functions is provided at the bottom of this notebook. We
highly recommend working in a `Jupyter notebook environment`_ in order to make use of
``alminer``'s visualization tools. We aim to keep adding new notebooks
relevant for various sub-fields in the future.

.. only:: latex

  .. note::
    To work with the tutorial notebook interactively |badge_pdf|

.. only:: html

  .. note::
    To work with the tutorial notebook interactively |badge_svg|


.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/1_query_tools
   tutorials/2_filter_explore
   tutorials/3_plot_results
   tutorials/4_create_reports
   tutorials/5_download_data
   tutorials/6_advanced_query

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: ALMA

   pages/scientific_categories
   pages/science_keywords

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: ALminer

   pages/query_keywords
   pages/api

Acknowledgements
----------------

``alminer`` has been developed through a collaboration between
`Allegro`_, the ALMA Regional Centre in
The Netherlands, and the University of Vienna as part of the `EMERGE-StG
project`_. This project has received
funding from the European Research Council (ERC) under the European
Union’s Horizon 2020 research and innovation programme (Grant agreement
No. 851435).

If you use ``alminer`` as part of your research, please consider citing
this `ASCL article`_ (ADS reference will
be added to the `GitHub`_ page when available).

``alminer`` makes use of different routines in
`Astropy`_ and `Astroquery`_. Please
also consider citing the following papers: - Astropy: `Astropy Collaboration et al. 2013`_
- Astroquery: `Ginsburg et al. 2019`_

We also acknowledge the work of Leiden University M.Sc. students, Robin
Mentel and David van Dop, who contributed to early versions of this
work.

Contact us
----------

If you encounter issues, please `open an issue`_.

If you have suggestions for improvement or would like to collaborate
with us on this project, please e-mail `Aida Ahmadi`_ and `Alvaro Hacar`_.

.. |badge_svg| image:: ../images/Binder_badge.svg
   :target: https://mybinder.org/v2/gh/emerge-erc/ALminer/main?filepath=notebooks/tutorial/ALminer_tutorial.ipynb

.. |badge_pdf| image:: ../images/Binder_badge.pdf
   :target: https://mybinder.org/v2/gh/emerge-erc/ALminer/main?filepath=notebooks/tutorial/ALminer_tutorial.ipynb

.. _GitHub: https://github.com/emerge-erc/ALminer
.. _open an issue: https://github.com/emerge-erc/ALminer/issues
.. _ASCL article: https://ascl.net/code/v/2971
.. _Aida Ahmadi: mailto:aahmadi@strw.leidenuniv.nl
.. _Alvaro Hacar: mailto:alvaro.hacar@univie.ac.at
.. _Ginsburg et al. 2019: https://ui.adsabs.harvard.edu/abs/2019AJ....157...98G/abstract
.. _Astropy Collaboration et al. 2013: https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A/abstract
.. _Astroquery: https://astroquery.readthedocs.io/en/latest/
.. _Astropy: https://www.astropy.org/
.. _EMERGE-StG project: https://emerge.alvarohacar.com
.. _Allegro: https://www.alma-allegro.nl/
.. _Jupyter notebook environment: https://jupyter.org/install
.. _tutorial Jupyter Notebook: https://nbviewer.jupyter.org/github/emerge-erc/ALminer/blob/main/notebooks/tutorial/ALminer_tutorial.ipynb?flush_cache=True
.. _ALMA science archive: https://almascience.eso.org/aq/
.. _pandas: https://pandas.pydata.org/
.. _pyvo: https://pyvo.readthedocs.io/en/latest