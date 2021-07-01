"""
ALminer
=======
ALMA archive mining and visualization toolkit
"""
##############################################
# Libraries
##############################################
import os
import re
import pandas as pd
import numpy as np
np.set_printoptions(threshold=np.inf)

# astropy
from astropy.coordinates import Angle
from astropy.coordinates import get_icrs_coordinates
from astropy.coordinates import name_resolve
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as const

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, NullFormatter
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# astroquery
from astroquery.alma import Alma

# TAP
from pyvo.dal import tap

##############################################
# Global Variables & Constants
##############################################

from .constants import band_names, band_color, band_min_freq, band_max_freq, \
    CO_line_names, CO_line_freq, CO_line_ha, CO_line_label, VALID_KEYWORDS_STR, \
    NEW_COLUMNS, COLUMN_TYPES


##############################################
# Setup Functions
##############################################

def _set_service():
    """Set the url for the ALMA TAP service."""
    service = tap.TAPService("https://almascience.eso.org/tap")
    return service


def _get_metadata():
    """Get TAP metadata in order to set units properly."""
    service = _set_service()
    metadata_query = "SELECT column_name, datatype, unit, ucd, utype, description from TAP_SCHEMA.columns"
    TAP_metadata = service.search(metadata_query)
    return pd.DataFrame(TAP_metadata).set_index('column_name')


##############################################
# Query Functions
##############################################

def conesearch(ra, dec, search_radius=1., point=True, public=True, published=None, print_targets=True,
               print_query=False):
    """
    Query the ALMA archive for a given position and radius around it.

    Parameters
    ----------
    ra : float
         Right ascension in degrees (ICRS).
    dec : float
         Declination in degrees (ICRS).
    search_radius : float, optional
         (Default value = 1. arcmin)
         Search radius (in arcmin) around the source coordinates.
    point : bool, optional
         (Default value = True)
         Search whether the phase center of the observations is contained within the search_radius (point=True)
         or whether any part of the observed region overlaps with the cone extending the search_radius (point=False).
         Note that point=True is much faster than point=False but the latter should be used if one is interested in
         searching for mosaics.
    public : bool, optional
         (Default value = True)
         Search for public data (public=True), proprietary data (public=False),
         or both public and proprietary data (public=None).
    published : bool, optional
         (Default value = None)
         Search for published data only (published=True), unpublished data only (published=False),
         or both published and unpublished data (published=None).
    print_query : bool, optional
         (Default value = True)
         Print the ADQL TAP query to the terminal.
    print_targets : bool, optional
         (Default value = False)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing the query results

    """
    search_radius = search_radius * u.arcmin
    if point:
        query = "SELECT * FROM ivoa.obscore WHERE 1 = CONTAINS(POINT('ICRS', s_ra, s_dec), " \
                "CIRCLE('ICRS',{},{},{}))".format(ra, dec, search_radius.to(u.deg).value)
    else:
        query = "SELECT * FROM ivoa.ObsCore WHERE 1 = INTERSECTS(CIRCLE('ICRS',{},{},{}), " \
                "s_region)".format(ra, dec, search_radius.to(u.deg).value)

    if public:
        query = "{} AND data_rights LIKE '%Public%'".format(query)
    elif not public and public is not None:
        query = "{} AND data_rights LIKE '%Proprietary%'".format(query)

    if print_query:
        print("Your query is: {}".format(query))

    TAP_df = run_query(query)
    if TAP_df is not None:
        if published:  # case pf published = True
            TAP_df = TAP_df[TAP_df['publication_year'].notnull()]
        elif not published and published is not None:  # case of published = False
            TAP_df = TAP_df[TAP_df['publication_year'].isnull()]
        filtered_df = filter_results(TAP_df, print_targets=print_targets)
        return filtered_df


def run_query(query_str):
    """
    Run the TAP query through PyVO service.

    Parameters
    ----------
    query_str : str
         ADQL query to send to the PyVO TAP service

    Returns
    -------
    pandas.DataFrame containing the query results

    """
    service = _set_service()
    # Run query
    pyvo_TAP_results = service.search(query_str)  # for large queries add maxrec=1000000
    # Transform output into astropy table first, then to a pandas DataFrame
    TAP_df = pyvo_TAP_results.to_table().to_pandas()
    # the column publication_year must be in 'object' type because it contains numbers and NaNs
    TAP_df['publication_year'] = TAP_df['publication_year'].astype('object')
    return TAP_df


def target(sources, search_radius=1., point=True, public=True, published=None, print_query=False, print_targets=True):
    """
    Query targets by name.

    This is done by using the astropy SESAME resolver to get the target's coordinates and then the ALMA archive
    is queried for those coordinates and a search_radius around them. The SESAME resolver searches multiple databases
    (Simbad, NED, VizieR) to parse names commonly found throughout literature and returns their coordinates. If the
    target is not resolved in any of these databases, consider using the 'keysearch' function and query the archive
    using the 'target_name' keyword (e.g. keysearch({'target_name': sources})).

    Parameters
    ----------
    sources : list of str
         list of sources by name.
         (IMPORTANT: source names must be identified by at least one of Simbad, NED, or Vizier)
    search_radius : float, optional
         (Default value = 1. arcmin)
         Search radius (in arcmin) around the source coordinates.
    point : bool, optional
         (Default value = True)
         Search whether the phase center of the observations is contained within the search_radius (point=True)
         or whether any part of the observed region overlaps with the cone extending the search_radius (point=False).
         Note that point=True is much faster than point=False but the latter should be used if one is interested in
         searching for mosaics.
    public : bool, optional
         (Default value = True)
         Search for public data (public=True), proprietary data (public=False),
         or both public and proprietary data (public=None).
    published : bool, optional
         (Default value = None)
         Search for published data only (published=True), unpublished data only (published=False),
         or both published and unpublished data (published=None).
    print_query : bool, optional
         (Default value = True)
         Print the ADQL TAP query to the terminal.
    print_targets : bool, optional
         (Default value = False)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing the query results.

    See Also
    --------
    keysearch : Query the ALMA archive for any (string-type) keywords defined in ALMA TAP system.

    """
    print("================================")
    print("alminer.target results ")
    print("================================")
    complete_results = []
    # go through list of sources provided by user and add query results to a list
    for s in sources:
        print("Target = {}".format(s))
        try:
            # Get source coodinates from astropy SESAME resolver querying multiple databases (SIMBAD, NED, Vizier)
            source_pos = get_icrs_coordinates(s)
            TAP_df = conesearch(ra=source_pos.ra.deg, dec=source_pos.dec.deg, search_radius=search_radius,
                                point=point, public=public, published=published, print_query=print_query,
                                print_targets=print_targets)
            if TAP_df is not None:
                complete_results.append(TAP_df)
        except name_resolve.NameResolveError as err:  # source coords not found in SESAME resolver
            print(err)
            print("Try keysearch function instead: keysearch({{'target_name':['{}']}}).".format(s))
            print("--------------------------------")
            pass
    # if the list of query results is not empty, concatenate them together into one DataFrame
    if complete_results:
        obs = pd.concat(complete_results)
        # need to reset the index of DataFrame so the indices in the final DataFrame are consecutive
        obs = obs.reset_index(drop=True)
        return obs
    else:
        print("No observations found for any sources in this list.")
        print("--------------------------------")


def catalog(target_df, search_radius=1., point=True, public=True, published=None, print_query=False,
            print_targets=True):
    """
    Query the ALMA archive for a list of coordinates or a catalog of sources based on their coordinates.

    Parameters
    ----------
    target_df : pandas.DataFrame
         Source names and coordinates.
         Index:
            RangeIndex
         Columns:
            Name: Name, dtype: str, description: target name (can be numbers or dummy names)
            Name: RAJ2000, dtype: float64, description: right ascension in degrees (ICRS)
            Name: DEJ2000, dtype: float64, description: declination in degrees (ICRS)
        
    search_radius : float, optional
         (Default value = 1. arcmin)
         Search radius (in arcmin) around the source coordinates.
    point : bool, optional
         (Default value = True)
         Search whether the phase center of the observations is contained within the search_radius (point=True)
         or whether any part of the observed region overlaps with the cone extending the search_radius (point=False).
         Note that point=True is much faster than point=False but the latter should be used if one is interested in
         searching for mosaics.
    public : bool, optional
         (Default value = True)
         Search for public data (public=True), proprietary data (public=False),
         or both public and proprietary data (public=None).
    published : bool, optional
         (Default value = None)
         Search for published data only (published=True), unpublished data only (published=False),
         or both published and unpublished data (published=None).
    print_query : bool, optional
         (Default value = True)
         Print the ADQL TAP query to the terminal.
    print_targets : bool, optional
         (Default value = False)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing the query results.

    """
    print("================================")
    print("alminer.catalog results")
    print("================================")
    complete_results = []
    for p in range(target_df.shape[0]):
        print("Target = {}".format(target_df.Name[p]))
        source_pos = SkyCoord(target_df.RAJ2000[p] * u.deg, target_df.DEJ2000[p] * u.deg, frame='icrs')
        TAP_df = conesearch(ra=source_pos.ra.deg, dec=source_pos.dec.deg, search_radius=search_radius,
                            point=point, public=public, published=published, print_query=print_query,
                            print_targets=print_targets)
        if TAP_df is not None:
            complete_results.append(TAP_df)
    # if the list of query results is not empty, concatenate them together into one DataFrame
    if complete_results:
        obs = pd.concat(complete_results)
        # need to reset the index of DataFrame so the indices in the final DataFrame are consecutive
        obs = obs.drop_duplicates().reset_index(drop=True)
        return obs
    else:
        print("No observations found for any sources in this catalog.")
        print("--------------------------------")


def keysearch(search_dict, public=True, published=None, print_query=False, print_targets=True):
    """
    Query the ALMA archive for any (string-type) keywords defined in ALMA TAP system.

    Parameters
    ----------
    search_dict : dict[str, list of str]
         Dictionary of keywords in the ALMA archive and their values. Values must be formatted as a list.
         A list of valid keywords are stored in VALID_KEYWORDS_STR variable.
    public : bool, optional
         (Default value = True)
         Search for public data (public=True), proprietary data (public=False),
         or both public and proprietary data (public=None).
    published : bool, optional
         (Default value = None)
         Search for published data only (published=True), unpublished data only (published=False),
         or both published and unpublished data (published=None).
    print_query : bool, optional
         (Default value = True)
         Print the ADQL TAP query to the terminal.
    print_targets : bool, optional
         (Default value = False)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing the query results.

    Notes
    -----
    The power of this function is in combining keywords. When multiple keywords are provided, they are
    queried using 'AND' logic, but when multiple values are provided for a given keyword, they are queried using
    'OR' logic. If a given value contains spaces, its constituents are queried using 'AND' logic. The only exception
    to this rule is the 'target_name' keyword.

    Examples
    --------
    keysearch({"proposal_abstract": ["high-mass star formation outflow disk"]})
         will query the archive for projects with the words
         "high-mass" AND "star" AND "formation" AND "outflow" AND "disk" in their proposal abstracts.

    keysearch({"proposal_abstract": ["high-mass", "star", "formation", "outflow", "disk"]})
         will query the archive for projects with the words
         "high-mass" OR "star" OR "formation" OR "outflow" OR "disk" in their proposal abstracts.

    keysearch({"proposal_abstract": ["star formation"], "scientific_category":['Galaxies']})
         will query the archive for projects with the words
         "star" AND "formation" in their proposal abstracts AND
         projects that are within the scientific_category of 'Galaxies'.

    """
    print("================================")
    print("alminer.keysearch results ")
    print("================================")
    # Add keyword to the query dictionary for the data rights (Public, Proprietary, or both)
    if public:
        search_dict['data_rights'] = ['Public']
    elif not public and public is not None:
        search_dict['data_rights'] = ['Proprietary']
    # Add scan intent keyword to the query dictionary to be the science target by default
    search_dict['scan_intent'] = ['TARGET']

    # compile a list of queries based on all keywords provided
    full_query_list = []
    for keyword, values in search_dict.items():
        # Catch if a wrong keyword is used and give appropriate error
        assert keyword in VALID_KEYWORDS_STR, "Invalid keyword, must be one of: {}".format(VALID_KEYWORDS_STR)
        # Convert underscores and spaces in the target name to wildcard
        # target_name is always queried with OR logic
        if keyword == 'target_name':
            values = [v.replace('_', '%') for v in values]
            values = [v.replace(' ', '%') for v in values]
            # Create queries for a given keyword using 'OR' logic between different values and accounting for
            # the case-sensitivity
            current_query = ["LOWER({}) LIKE '%{}%'".format(keyword, v.lower()) for v in values]
            full_query_list.append("({})".format(" OR ".join(current_query)))
        # Account for AND/OR logic for keywords that are not target_name
        else:
            keyword_query_list = []
            for v in values:
                # If there are spaces in the values of a given keyword, split them out and query them with AND logic
                if re.search(r"\s", v):
                    split_values = v.split()
                    current_query = ["LOWER({}) LIKE '%{}%'".format(keyword, s.lower()) for s in split_values]
                    keyword_query_list.append("({})".format(" AND ".join(current_query)))
                # If separate words are provided as values, query them with OR logic
                else:
                    keyword_query_list.append("LOWER({}) LIKE '%{}%'".format(keyword, v.lower()))
            full_query_list.append("({})".format(" OR ".join(keyword_query_list)))
    # Put together the entire query with 'AND' logic between different keywords
    full_query = "SELECT * FROM ivoa.obscore WHERE {} ORDER BY proposal_id".format(" AND ".join(full_query_list))
    if print_query:
        print("Your query is: {}".format(full_query))
    TAP_df = run_query(full_query)
    # Filter whether the user wants published data, unpublished data, or both (default)
    if published:  # case pf published = True
        TAP_df = TAP_df[TAP_df['publication_year'].notnull()]
    elif not published and published is not None:  # case pf published = False
        TAP_df = TAP_df[TAP_df['publication_year'].isnull()]
    return filter_results(TAP_df, print_targets=print_targets)


def _get_freq_res(frequency_support, freq_min, freq_max, em_res_power):
    """
    Given the minimum and maximum frequency, parse the 'frequency_support' and extract the frequency resolution.
    """
    for s in frequency_support.split('[')[1:]:
        if s.startswith("{:.2f}..{:.2f}".format(freq_min, freq_max)):
            return float(s.split(',')[1].replace('kHz', ''))
    # if the parsing is not successful, calculate the frequency resolution manually
    return float((freq_max * u.GHz / em_res_power).to(u.kHz).value.round(decimals=2))


def _set_plot_params(ax, title, xlabel, ylabel, legend=True, log=False):
    """Set the plot parameters like titles, labels, and scaling."""
    ax.autoscale(enable=True, axis='both', tight=False)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if log:
        ax.set_yscale('log')
        ax.set_ylim(ymin=1e-1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax.yaxis.set_minor_formatter(NullFormatter())
    if legend:
        ax.legend(fontsize=15, loc='best')
        ax.autoscale(enable=True, axis='both', tight=False)
        ax.set_ylim(ymin=1e-1)


def filter_results(TAP_df, print_targets=True):
    """
    Add a few new useful columns to the pandas.DataFrame with the query results from the PyVO TAP service and
    return the full query DataFrame and optionally a summary of the results.

    Parameters
    ----------
    TAP_df : pandas.DataFrame
         This is likely the output of 'run_query' function.
    print_targets : bool, optional
         (Default value = True)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing the query results.

    """
    # add new columns to the DataFrame
    TAP_COLUMNS = TAP_df.columns.tolist()
    data = TAP_df.reindex(NEW_COLUMNS + TAP_COLUMNS, axis=1)
    data = data.astype(COLUMN_TYPES)
    data = data.reset_index(drop=True)  # needed to renumber the index of the rows for looping over below
    if not data.empty:
        # calculate the relevant variables
        data['min_freq_GHz'] = [(em_max * u.m).to(u.GHz, equivalencies=u.spectral()).value for em_max in data['em_max']]
        data['max_freq_GHz'] = [(em_min * u.m).to(u.GHz, equivalencies=u.spectral()).value for em_min in data['em_min']]
        data['central_freq_GHz'] = data[['min_freq_GHz', 'max_freq_GHz']].mean(axis=1)
        data['bandwidth_GHz'] = data['max_freq_GHz'] - data['min_freq_GHz']
        data['line_sens_10kms'] = data['sensitivity_10kms'] * 1000.0  # in uJy/beam
        data['cont_sens_bandwidth'] = data['cont_sensitivity_bandwidth'] * 1000.0  # in uJy/beam

        # parse the 'frequency_support' string to determine the frequency_resolution EXACTLY
        data['freq_res_kHz'] = [_get_freq_res(fs, data['min_freq_GHz'][i], data['max_freq_GHz'][i],
                                              data['em_res_power'][i])
                                for i, fs in enumerate(data['frequency_support'])]

        nchan = data['bandwidth_GHz'] * 1000000.0 / data['freq_res_kHz']
        data['vel_res_kms'] = const.c.to(u.km / u.s).value / data['em_res_power']
        data['line_sens_native'] = data['line_sens_10kms'] * np.sqrt(10.0 / data['vel_res_kms'] / nchan)

        # add values to the table in desired precision
        data['Obs'] = [i+1 for i in np.arange(data.shape[0])]
        data['project_code'] = data['proposal_id']
        data['ALMA_source_name'] = data['target_name']
        data['RAJ2000'] = data['s_ra']
        data['DEJ2000'] = data['s_dec']
        data['ang_res_arcsec'] = data['s_resolution'].round(decimals=3)
        data['min_freq_GHz'] = data['min_freq_GHz'].round(decimals=2)
        data['max_freq_GHz'] = data['max_freq_GHz'].round(decimals=2)
        data['central_freq_GHz'] = data['central_freq_GHz'].round(decimals=2)
        data['bandwidth_GHz'] = data['bandwidth_GHz'].round(decimals=3)
        data['freq_res_kHz'] = data['freq_res_kHz'].round(decimals=2)
        data['vel_res_kms'] = data['vel_res_kms'].round(decimals=3)
        data['LAS_arcsec'] = data['spatial_scale_max'].round(decimals=3)
        data['FoV_arcsec'] = (data['s_fov'] * 3600.).round(decimals=3)
        data['cont_sens_bandwidth'] = data['cont_sens_bandwidth'].round(decimals=3)
        data['line_sens_10kms'] = data['line_sens_10kms'].round(decimals=2)
        data['line_sens_native'] = data['line_sens_native'].round(decimals=2)
        data['MOUS_id'] = data['member_ous_uid']
        data = data.drop_duplicates().reset_index(drop=True)
        summary(data, print_targets=print_targets)
        return data
    else:
        print("--------------------------------")
        print("No observations found.")
        print("--------------------------------")


##############################################
# Analysis Functions
##############################################

def summary(observations, print_targets=True):
    """
    Print a summary of the observations.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    print_targets : bool, optional
         (Default value = True)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    """
    if not observations.empty:
        print("--------------------------------")
        unique_projects = observations.proposal_id.unique()
        print("Number of projects = {}".format(len(unique_projects)))
        unique_observations = observations.drop_duplicates(subset=['group_ous_uid', 'member_ous_uid', 'target_name'],
                                                           ignore_index=True)
        print("Number of observations = {}".format(np.shape(unique_observations)[0]))
        unique_subbands = observations.drop_duplicates(subset=['min_freq_GHz', 'max_freq_GHz', 'freq_res_kHz'],
                                                       ignore_index=True)
        print("Number of unique subbands = {}".format(np.shape(unique_subbands)[0]))
        print("Total number of subbands = {}".format(np.shape(observations)[0]))
        # Filter out the repeating source names
        observed_sources = observations.target_name.unique().tolist()
        if print_targets:
            print("{} target(s) with ALMA data = {}".format(len(observed_sources), observed_sources))
        else:
            print("Total number of targets with ALMA data = {}".format(len(observed_sources)))
        print("--------------------------------")
    else:
        print("--------------------------------")
        print("No observations found.")
        print("--------------------------------")


def explore(observations, allcols=False, allrows=False):
    """
    Control how much of the pandas.DataFrame with the query results is presented in the displayed table.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    allcols : bool, optional
         (Default value = False)
         Show all 81 columns (allcols=True), or the first 18 columns (allcols=False).
    allrows : bool, optional
         (Default value = False)
         Show all rows in the DataFrame (allrows=True), or just a summary (allrows=False).

    Returns
    -------
    pandas.DataFrame containing the query results displayed to the user interface as specified by the user.

    """
    pd.set_option('display.max_rows', 0)
    pd.set_option('display.max_columns', 0)
    if allrows:
        pd.set_option('display.max_rows', None)
    if not allcols:
        observations = observations.iloc[:, : len(NEW_COLUMNS)]
    return observations
    

def get_units(column):
    """Print the units for a given column in the query results DataFrame."""
    COLUMN_UNITS = {'Obs': "''", 'project_code': "''", 'ALMA_source_name': "''", 'RAJ2000': 'deg', 'DEJ2000': 'deg',
                    'ang_res_arcsec': 'arcsec', 'min_freq_GHz': 'GHz', 'max_freq_GHz': 'GHz', 'central_freq_GHz': 'GHz',
                    'bandwidth_GHz': 'GHz', 'freq_res_kHz': 'kHz', 'vel_res_kms': 'km/s', 'LAS_arcsec': 'arcsec',
                    'FoV_arcsec': 'arcsec', 'cont_sens_bandwidth': 'uJy/beam', 'line_sens_10kms': 'uJy/beam',
                    'line_sens_native': 'uJy/beam', 'MOUS_id': "''"}
    try:
        TAP_metadata = _get_metadata()
        return TAP_metadata.loc[column, 'unit']
    except KeyError:
        return COLUMN_UNITS.get(column, "No units found.")


def get_description(column):
    """Print the description of a given column in the query results DataFrame."""
    TAP_metadata = _get_metadata()
    COLUMN_DESCRIPTION = {'Obs': "Observation number in the DataFrame",
                          'project_code': TAP_metadata.loc['proposal_id', 'description'],
                          'ALMA_source_name': TAP_metadata.loc['target_name', 'description'],
                          'RAJ2000': TAP_metadata.loc['s_ra', 'description'],
                          'DEJ2000': TAP_metadata.loc['s_dec', 'description'],
                          'ang_res_arcsec': TAP_metadata.loc['s_resolution', 'description'],
                          'min_freq_GHz': "Minimum frequency of the subband in GHz",
                          'max_freq_GHz': "Maximum frequency of the subband in GHz",
                          'central_freq_GHz': "Central frequency of the subband in GHz",
                          'bandwidth_GHz': "Bandwidth of the subband in GHz",
                          'freq_res_kHz': "Spectral resolution of the subband in kHz",
                          'vel_res_kms': "Velocity resolution of the subband in km/s",
                          'LAS_arcsec': "Largest recoverable angular scales in arcsec",
                          'FoV_arcsec': TAP_metadata.loc['s_fov', 'description'],
                          'cont_sens_bandwidth': TAP_metadata.loc['cont_sensitivity_bandwidth', 'description'],
                          'line_sens_10kms': TAP_metadata.loc['sensitivity_10kms', 'description'],
                          'line_sens_native': "Estimated noise in subband's native spectral resolution",
                          'MOUS_id': TAP_metadata.loc['member_ous_uid', 'description']}
    try:
        return TAP_metadata.loc[column, 'description']
    except KeyError:
        return COLUMN_DESCRIPTION.get(column, "No description found.")


def get_info(column):
    """Print the description and units of a given column in the query results DataFrame.

    Parameters
    ----------
    column : str
         A column in the pandas.DataFrame query table.

    """
    print("--------------------------------")
    print("Column: {}".format(column))
    print("--------------------------------")
    print("Description: {}".format(get_description(column)))
    print("Units: {}".format(get_units(column)))
    print("--------------------------------")


def line_coverage(observations, line_freq, z=0., line_name='', print_summary=True, print_targets=True):
    """
    Determine how many observations were observed at a given frequency (+redshift).

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    line_freq : float64
         Frequency of the line of interest in GHz.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequency given in 'line_freq' parameter should be shifted.
    line_name : str, optional
         (Default value = '')
         Name of the line specified in 'line_freq'.
    print_summary : bool, optional
         (Default value = True)
         Print a summary of the observations to the terminal.
    print_targets : bool, optional
         (Default value = True)
         Print a list of targets with ALMA data (ALMA source names) to the terminal.

    Returns
    -------
    pandas.DataFrame containing all observations of line of interest.

    """
    # shift the line by redshift
    line_freq_shifted = line_freq * (1. / (1. + z))
    covered_df = observations[(observations["min_freq_GHz"] < line_freq_shifted)
                              & (observations["max_freq_GHz"] > line_freq_shifted)]
    if print_summary:
        print("--------------------------------")
        if z:
            print("Summary of '{}' observations at {} GHz ({} GHz at z={})".format(line_name, line_freq,
                                                                                   round(line_freq_shifted, 3), z))
        else:
            print("Summary of '{}' observations at {} GHz".format(line_name, round(line_freq_shifted, 3)))
        summary(covered_df, print_targets=print_targets)
    return covered_df


def CO_lines(observations, z=0., print_summary=True, print_targets=True):
    """
    Determine how many CO, 13CO, and C18O lines were observed in the provided query DataFrame.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequencies should be shifted.
    print_summary : bool, optional
         (Default value = True)
         Print a summary of the observations for each (redshifted) CO, 13CO, and C18O line to the terminal.
    print_targets : bool, optional
         (Default value = True)
         Print the target names (ALMA source names) with ALMA data for each (redshifted) CO, 13CO, and C18O line
         to the terminal.

    Returns
    -------
    pandas.DataFrame containing all observations of (redshifted) CO, 13CO, and C18O lines.

    """
    CO_df_list = []
    for t, line in enumerate(CO_line_names):
        line_df = line_coverage(observations, line_freq=CO_line_freq[line], z=z,
                                line_name=CO_line_names[t], print_summary=print_summary, print_targets=print_targets)
        if not line_df.empty:
            CO_df_list.append(line_df)
    if CO_df_list:
        CO_df = pd.concat(CO_df_list)
        # need to reset the index of DataFrame so the indices in the final DataFrame are consecutive
        CO_df = CO_df.drop_duplicates().reset_index(drop=True)
        return CO_df
    else:
        print("Found no ALMA observations covering transitions of CO, 13CO, or C18O.")
        print("--------------------------------")


##############################################
# Plotting Functions
##############################################

def _draw_freq(ax, direction, mark_freq, z=0.):
    """Draw a horizontal/vertical dashed line at a given (redshifted) frequency on the requested axis."""
    # if the value is not given as a list, turn it into a list and loop over it to draw lines
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if not isinstance(mark_freq, list):
        mark_freq = [mark_freq]
    if isinstance(mark_freq, list):
        # convert items in mark_freq to float type in case they are not already
        mark_freq = [float(item) for item in mark_freq]
        for f, freq in enumerate(mark_freq):
            freq_shifted = freq * (1. / (1. + z))
            if direction == 'vertical':
                if xmin <= freq_shifted <= xmax:
                    ax.axvline(freq_shifted, linestyle="--", color="black", zorder=3)
            elif direction == 'horizontal':
                if ymin <= freq_shifted <= ymax:
                    ax.axhline(freq_shifted, linestyle="--", color="black", zorder=3)
    else:
        print("Please provide a valid list of frequencies in GHz. For example: mark_freq=[220.5, 345.0]")


def _draw_CO(ax, direction, z=0.):
    """
    Draw a horizontal/vertical dashed line at (redshifted) frequencies of CO, 13CO, and C18O lines on the requested
    axis.

    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    corners_list = []
    for t, line in enumerate(CO_line_names):
        freq = CO_line_freq[line]
        freq_shifted = freq * (1. / (1. + z))
        if direction == 'vertical':
            if xmin <= freq_shifted <= xmax:
                ax.vlines(freq_shifted, 1e-1, 3*ymax, linestyle="dotted", color="black", zorder=1)

                annotation = ax.annotate(CO_line_label[line], xy=(freq_shifted, 3.5*ymax),
                                         xycoords='data', rotation=direction, ha=CO_line_ha[line],
                                         fontsize=10, zorder=3, clip_on=True)
                ax.set_ylim(1e-1, 50*ymax)
                bbox = annotation.get_window_extent(plt.gcf().canvas.get_renderer())
                bbox_data = bbox.transformed(ax.transData.inverted())
                corners_list.append(bbox_data.corners()[:, 1])
        elif direction == 'horizontal':
            if ymin <= freq_shifted <= ymax:
                if CO_line_ha[line] == 'left':
                    line_va = 'bottom'
                elif CO_line_ha[line] == 'right':
                    line_va = 'top'
                else:
                    line_va = 'center'
                ax.hlines(freq_shifted, xmin, xmax, linestyle="dotted", color="black", zorder=3)
                annotation = ax.text(xmax, freq_shifted, CO_line_label[line], rotation=direction,
                                     va=line_va, zorder=3, fontsize=10)
                bbox = annotation.get_window_extent(plt.gcf().canvas.get_renderer())
                bbox_data = bbox.transformed(ax.transData.inverted())
                ax.update_datalim(bbox_data.corners())
                ax.set_xlim(xmin=xmin)
    if len(corners_list) > 0:
        ax.set_ylim(ymax=np.max(corners_list))


def _band_detection(observations, min_freq, max_freq):
    """
    Obtain a subset of the provided pandas.DataFrame that contains observations in the range of minimum frequency
    (min_freq) to maximum frequency (max_freq) in GHz. This is used to define the frequency range of various ALMA bands.

    """
    band_obs = observations[(observations['central_freq_GHz'] >= min_freq)
                            & (observations['central_freq_GHz'] <= max_freq)]
    return band_obs


def plot_overview(observations, mark_freq='', z=0., mark_CO=False, showfig=True, savefig=None):
    """
    Create overview plots of observed frequencies, angular resolution, LAS, frequency and velocity resolutions.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    mark_freq : list of float64, optional
         (Default value = '')
         A list of frequencies to mark on the plot with dashed lines.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequencies given in 'mark_freq' and 'mark_CO' parameters should be shifted.
         Currently only one redshift can be given for all targets.
    mark_CO : bool, optional
         (Default value = False)
         Mark CO, 13CO, and C18O frequencies on the plot with dashed lines.
    showfig : bool, optional
         (Default value = True)
         Display the plot (showfig=True) or not (showfig=False).
    savefig : str, optional
         (Default value = None)
         Filename (without an extension) for the plot to be saved as.
         Default file extension is PDF. Figure is saved in a subdirectory called 'reports' within the current
         working directory. If the directory doesn't exist, it will be created. Default quality is dpi=300.

    """
    if not observations.empty:
        # Overview plots of observed frequencies, coloured per band
        fig = plt.figure(figsize=(12, 15))
        ax0 = plt.subplot(3, 2, (1, 2))
        bin_width = 5
        for b, band in enumerate(band_names):
            # split observations into different bands
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                # Set the binning to be done in 5 GHz widths for the overview plot
                bins = np.arange(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), bin_width)
                # If the binning is too large, then set the binning based on the data
                if len(bins) < 5:
                    bins = np.linspace(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), 20)
                ax0.hist(band_df['central_freq_GHz'], bins=bins, label="{}: obs. = {}"
                         .format(band, band_df.shape[0]),
                         color=band_color[band], edgecolor='black')
        _set_plot_params(ax=ax0, title="Observed Frequencies", xlabel="Frequency (GHz)",
                         ylabel="Number of observations", log=True)
        # draw vertical lines at a given frequency as requested by user
        if mark_freq:
            _draw_freq(ax=ax0, direction='vertical', mark_freq=mark_freq, z=z)
        # draw vertical lines at frequencies associated with CO, 13CO, & C18O
        if mark_CO:
            _draw_CO(ax=ax0, direction='vertical', z=z)

        # Resolution
        ax1 = plt.subplot(3, 2, 3)
        bins = np.linspace(min(observations['ang_res_arcsec']), max(observations['ang_res_arcsec']), 20)
        ax1.hist(observations['ang_res_arcsec'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['ang_res_arcsec']), 2), round(max(observations['ang_res_arcsec']), 2)),
                 color="lightgrey", edgecolor='black')
        _set_plot_params(ax=ax1, title="Angular resolution", xlabel="Beam size (arcsec)",
                         ylabel="Number of observations", log=True)

        # Largest Angular Scale
        ax2 = plt.subplot(3, 2, 4)
        bins = np.linspace(min(observations["LAS_arcsec"]), max(observations["LAS_arcsec"]), 20)
        ax2.hist(observations['LAS_arcsec'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['LAS_arcsec']), 2), round(max(observations['LAS_arcsec']), 2)),
                 color="lightgrey", edgecolor='black')
        _set_plot_params(ax=ax2, title="Largest Angular Scale", xlabel="LAS (arcsec)",
                         ylabel="Number of observations", log=True)

        # Freq. Resolution
        ax3 = plt.subplot(3, 2, 5)
        bins = np.linspace(min(observations['freq_res_kHz']), max(observations['freq_res_kHz']), 20)
        ax3.hist(observations['freq_res_kHz'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['freq_res_kHz']), 2), round(max(observations['freq_res_kHz']), 2)),
                 color="lightgrey", edgecolor='black')
        _set_plot_params(ax=ax3, title="Frequency resolution", xlabel="Frequency resolution (kHz)",
                         ylabel="Number of observations", log=True)

        # Velocity Resolution
        ax4 = plt.subplot(3, 2, 6)
        bins = np.linspace(min(observations['vel_res_kms']), max(observations['vel_res_kms']), 20)
        ax4.hist(observations['vel_res_kms'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['vel_res_kms']), 2), round(max(observations['vel_res_kms']), 2)),
                 color="lightgrey", edgecolor='black')
        _set_plot_params(ax=ax4, title="Velocity resolution", xlabel="Velocity resolution (km s$^{-1}$)",
                         ylabel="Number of observations", log=True)

        fig.tight_layout()
        plt.rcParams['figure.dpi'] = 300
        _save_plot(savefig)
        if showfig:
            plt.show()
        plt.close()
    else:
        print("--------------------------------")
        print("Nothing to plot: the DataFrame provided is empty.")
        print("--------------------------------")


def plot_bands(observations, mark_freq='', z=0., mark_CO=False, showfig=True, savefig=None):
    """
    Create overview and detailed plots of observed frequencies in each band.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    mark_freq : list of float64, optional
         (Default value = '')
         A list of frequencies to mark on the plot with dashed lines.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequencies given in 'mark_freq' and 'mark_CO' parameters should be shifted.
         Currently only one redshift can be given for all targets.
    mark_CO : bool, optional
         (Default value = False)
         Mark CO, 13CO, and C18O frequencies on the plot with dashed lines.
    showfig : bool, optional
         (Default value = True)
         Display the plot (showfig=True) or not (showfig=False).
    savefig : str, optional
         (Default value = None)
         Filename (without an extension) for the plot to be saved as.
         Default file extension is PDF. Figure is saved in a subdirectory called 'reports' within the current
         working directory. If the directory doesn't exist, it will be created. Default quality is dpi=300.

    """
    if not observations.empty:
        # Overview plots of observed frequencies, coloured per band
        plt.figure(figsize=(12, 18))
        ax0 = plt.subplot(4, 2, (1, 2))
        bin_width = 5
        for b, band in enumerate(band_names):
            # split observations into different bands
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                # Set the binning to be done in 5 GHz widths for the overview plot
                bins = np.arange(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), bin_width)
                # If the binning is too large, then set the binning based on the data
                if len(bins) < 5:
                    bins = np.linspace(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), 20)
                ax0.hist(band_df['central_freq_GHz'], bins=bins, label="{}: obs. = {}"
                         .format(band, band_df.shape[0]), color=band_color[band], edgecolor='black', zorder=2)
        _set_plot_params(ax=ax0, title="Observed Frequencies", xlabel="Frequency (GHz)",
                         ylabel="Number of observations", log=True)
        # draw vertical lines at frequencies associated with CO, 13CO, & C18O
        if mark_CO:
            _draw_CO(ax=ax0, direction='vertical', z=z)
        # draw vertical lines at a given frequency as requested by user
        if mark_freq:
            _draw_freq(ax=ax0, direction='vertical', mark_freq=mark_freq, z=z)

        # Zoom plots per band
        counter = 0
        for b, band in enumerate(band_names):
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                counter = counter + 1
                ax = plt.subplot(4, 2, 2 + counter)
                # Set the binning to be adjusted based on the data to 20 bins linearly spaced between the min/max
                bins = np.linspace(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), 20)
                ax.hist(band_df['central_freq_GHz'], bins=bins, color=band_color[band], edgecolor='black', zorder=2)
                _set_plot_params(ax=ax, title="{}: obs. = {}".format(band, band_df.shape[0]), xlabel="Frequency (GHz)",
                                 ylabel="Number of observations", log=True, legend=False)
                # draw vertical lines at a given frequency as requested by user
                if mark_freq:
                    _draw_freq(ax=ax, direction='vertical', mark_freq=mark_freq, z=z)
                # draw vertical lines at frequencies associated with CO, 13CO, & C18O
                if mark_CO:
                    _draw_CO(ax=ax, direction='vertical', z=z)

        plt.tight_layout()
        plt.rcParams['figure.dpi'] = 300
        _save_plot(savefig)
        if showfig:
            plt.show()
        plt.close()
    else:
        print("--------------------------------")
        print("Nothing to plot: the DataFrame provided is empty.")
        print("--------------------------------")
    

def plot_observations(observations, mark_freq='', z=0., mark_CO=False, showfig=True, savefig=None):
    """
    Create detailed plots of observations in each band. The x-axis displays the observation number 'Obs' column in
    the input DataFrame.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    mark_freq : list of float64, optional
         (Default value = '')
         A list of frequencies to mark on the plot with dashed lines.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequencies given in 'mark_freq' and 'mark_CO' parameters should be shifted.
         Currently only one redshift can be given for all targets.
    mark_CO : bool, optional
         (Default value = False)
         Mark CO, 13CO, and C18O frequencies on the plot with dashed lines.
    showfig : bool, optional
         (Default value = True)
         Display the plot (showfig=True) or not (showfig=False).
    savefig : str, optional
         (Default value = None)
         Filename (without an extension) for the plot to be saved as.
         Default file extension is PDF. Figure is saved in a subdirectory called 'reports' within the current
         working directory. If the directory doesn't exist, it will be created. Default quality is dpi=300.

    """
    if not observations.empty:
        # Find how many bands are needed to be plotted to set the number of rows in subplot
        nbands = 0
        for b, band in enumerate(band_names):
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                nbands = nbands + 1

        # Set the vertical size of plot based on number of bands detected
        plt.figure(figsize=(12, 5*nbands))
        counter = 0
        for b, band in enumerate(band_names):
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                counter = counter + 1
                ax = plt.subplot(nbands, 1, counter)
                ax.plot((band_df["Obs"], band_df["Obs"]), (band_df["min_freq_GHz"], band_df["max_freq_GHz"]),
                        color=band_color[band], linewidth=5.)
                _set_plot_params(ax=ax, title=band, xlabel="Observation Number", ylabel="Frequency (GHz)",
                                 legend=False, log=False)
                # draw horizontal lines at a given frequency as requested by user
                if mark_freq:
                    _draw_freq(ax=ax, direction='horizontal', mark_freq=mark_freq, z=z)
                # draw vertical lines at frequencies associated with CO, 13CO, & C18O
                if mark_CO:
                    _draw_CO(ax=ax, direction='horizontal', z=z)

        plt.tight_layout()
        plt.rcParams['figure.dpi'] = 300
        _save_plot(savefig)
        if showfig:
            plt.show()
        plt.close()
    else:
        print("--------------------------------")
        print("Nothing to plot: the DataFrame provided is empty.")
        print("--------------------------------")


def plot_line_overview(observations, line_freq, z=0., line_name='', showfig=True, savefig=None):
    """
    Create overview plots of observed frequencies, angular resolution, LAS, frequency and velocity resolutions,
    highlighting the observations of a give (redshifted) frequency with hatches on the bar plots.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    line_freq : float64
         Frequency of the line of interest in GHz.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequency given in 'line_freq' parameter should be shifted.
    line_name : str, optional
         (Default value = '')
         Name of the line specified in 'line_freq'.
    showfig : bool, optional
         (Default value = True)
         Display the plot (showfig=True) or not (showfig=False).
    savefig : str, optional
         (Default value = None)
         Filename (without an extension) for the plot to be saved as.
         Default file extension is PDF. Figure is saved in a subdirectory called 'reports' within the current
         working directory. If the directory doesn't exist, it will be created. Default quality is dpi=300.

    """
    if not observations.empty:
        # Overview plots of observed frequencies, coloured per band
        plt.figure(figsize=(12, 15))
        ax0 = plt.subplot(3, 1, 1)
        bin_width = 5
        # Get a DataFrame of the observations that include the requested line
        line_df = line_coverage(observations, line_freq=line_freq, z=z, line_name=line_name, print_summary=False,
                                print_targets=False)
        for b, band in enumerate(band_names):
            band_df = _band_detection(observations, min_freq=band_min_freq[band], max_freq=band_max_freq[band])
            if not band_df.empty:
                # Set the binning to be done in 5 GHz widths for the overview plot
                bins = np.arange(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), bin_width)
                # If the binning is too large, then set the binning based on the data
                if len(bins) < 5:
                    bins = np.linspace(min(band_df['central_freq_GHz']), max(band_df['central_freq_GHz']), 20)
                ax0.hist(band_df['central_freq_GHz'], bins=bins, label="{}: obs. = {}"
                         .format(band, band_df.shape[0]), color=band_color[band], edgecolor='black')
                if not line_df.empty:
                    if any(line_df.central_freq_GHz.isin(band_df.central_freq_GHz)):
                        ax0.hist(line_df['central_freq_GHz'], bins=bins, label="{} obs. = {}"
                                 .format(line_name, np.shape(line_df)[0]), color=band_color[band],
                                 hatch='///', edgecolor='black')
        _set_plot_params(ax=ax0, title="Observed Frequencies", xlabel="Frequency (GHz)",
                         ylabel="Number of observations", log=True)

        # Resolution
        ax1 = plt.subplot(3, 2, 3)
        bins = np.linspace(min(observations['ang_res_arcsec']), max(observations['ang_res_arcsec']), 20)
        ax1.hist(observations['ang_res_arcsec'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['ang_res_arcsec']), 2), round(max(observations['ang_res_arcsec']), 2)),
                 color="lightgrey", edgecolor='black')
        if not line_df.empty:
            ax1.hist(line_df['ang_res_arcsec'], bins=bins, label="min = {}, max = {}"
                     .format(round(min(line_df['ang_res_arcsec']), 2), round(max(line_df['ang_res_arcsec']), 2)),
                     color="lightgrey", hatch='///', edgecolor='black')
        _set_plot_params(ax=ax1, title="Angular resolution", xlabel="Beam size (arcsec)",
                         ylabel="Number of observations", log=True)

        # Largest Angular Scale
        ax2 = plt.subplot(3, 2, 4)
        bins = np.linspace(min(observations["LAS_arcsec"]), max(observations["LAS_arcsec"]), 20)
        ax2.hist(observations['LAS_arcsec'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['LAS_arcsec']), 2), round(max(observations['LAS_arcsec']), 2)),
                 color="lightgrey", edgecolor='black')
        if not line_df.empty:
            ax2.hist(line_df['LAS_arcsec'], bins=bins, label="min = {}, max = {}"
                     .format(round(min(line_df['LAS_arcsec']), 2), round(max(line_df['LAS_arcsec']), 2)),
                     color="lightgrey", hatch='///', edgecolor='black')
        _set_plot_params(ax=ax2, title="Largest Angular Scale", xlabel="LAS (arcsec)",
                         ylabel="Number of observations", log=True)

        # Freq. Resolution
        ax3 = plt.subplot(3, 2, 5)
        bins = np.linspace(min(observations['freq_res_kHz']), max(observations['freq_res_kHz']), 20)
        ax3.hist(observations['freq_res_kHz'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['freq_res_kHz']), 2), round(max(observations['freq_res_kHz']), 2)),
                 color="lightgrey", edgecolor='black')
        if not line_df.empty:
            ax3.hist(line_df['freq_res_kHz'], bins=bins, label="min = {}, max = {}"
                     .format(round(min(line_df['freq_res_kHz']), 2), round(max(line_df['freq_res_kHz']), 2)),
                     color="lightgrey", hatch='///', edgecolor='black')
        _set_plot_params(ax=ax3, title="Frequency resolution", xlabel="Frequency resolution (kHz)",
                         ylabel="Number of observations", log=True)

        # Velocity Resolution
        ax4 = plt.subplot(3, 2, 6)
        bins = np.linspace(min(observations['vel_res_kms']), max(observations['vel_res_kms']), 20)
        ax4.hist(observations['vel_res_kms'], bins=bins, label="min = {}, max = {}"
                 .format(round(min(observations['vel_res_kms']), 2), round(max(observations['vel_res_kms']), 2)),
                 color="lightgrey", edgecolor='black')
        if not line_df.empty:
            ax4.hist(line_df['vel_res_kms'], bins=bins, label="min = {}, max = {}"
                     .format(round(min(line_df['vel_res_kms']), 2), round(max(line_df['vel_res_kms']), 2)),
                     color="lightgrey", hatch='///', edgecolor='black')
        _set_plot_params(ax=ax4, title="Velocity resolution", xlabel="Velocity resolution (km s$^{-1}$)",
                         ylabel="Number of observations", log=True)

        plt.tight_layout()
        plt.rcParams['figure.dpi'] = 300
        _save_plot(savefig)
        if showfig:
            plt.show()
        plt.close()
    else:
        print("--------------------------------")
        print("Nothing to plot: the DataFrame provided is empty.")
        print("--------------------------------")


def plot_sky(observations, showfig=True, savefig=None):
    """
    Plot the distribution of the targets on the sky.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    showfig : bool, optional
         (Default value = True)
         Display the plot (showfig=True) or not (showfig=False).
    savefig : str, optional
         (Default value = None)
         Filename (without an extension) for the plot to be saved as.
         Default file extension is PDF. Figure is saved in a subdirectory called 'reports' within the current
         working directory. If the directory doesn't exist, it will be created. Default quality is dpi=300.

    """
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(111, projection="mollweide")
    ra = Angle(observations["RAJ2000"], unit=u.deg)
    ra = ra.wrap_at(180 * u.degree)
    dec = Angle(observations["DEJ2000"], unit=u.deg)
    ax.grid(True, linewidth=1, zorder=1)
    ax.set_xticklabels(['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h'])
    ax.scatter(ra.radian, dec.radian, color="#E5323B", zorder=2)
    ax.set_title("Sky distribution", fontsize=18)
    plt.rcParams['figure.dpi'] = 300
    _save_plot(savefig)
    if showfig:
        plt.show()
    plt.close()


##############################################
# Writing Files & Plots
##############################################

def _save_plot(savefig):
    """Save the current plot in the 'reports' subdirectory in the current working directory.
    If the directory doesn't exist, it will be created.

    """
    # Default path and file suffix for saving plots
    reports_path = './reports'
    suffix = '.pdf'
    if savefig is not None:
        # Create the reports directory if it does not exist
        if not os.path.isdir(reports_path):
            os.makedirs(reports_path)
        plt.savefig(os.path.join(reports_path, savefig + suffix), dpi=300, bbox_inches='tight')


def save_table(observations, filename="mytable"):
    """Write the DataFrame with the query results to a table in CSV format.

    The table will be saved in the 'tables' subdirectory within the current working directory.
    If the directory doesn't exist, it will be created.

    Parameters
    ----------
    observations : pandas.DataFrame
    filename : str
         (Default value = "mytable")
         Name of the table to be saved in the 'tables' subdirectory.

    """
    tables_path = './tables'
    suffix = '.csv'
    if not os.path.isdir(tables_path):
        os.makedirs(tables_path)
    observations.to_csv(os.path.join(tables_path, filename + suffix), index=False)


def save_source_reports(observations, mark_freq='', z=0., mark_CO=False):
    """
    Create overview plots of observed frequencies, angular resolution, LAS, frequency and velocity resolutions for
    each source in the provided DataFrame and save them in PDF format in the 'reports' subdirectory.
    If the directory doesn't exist, it will be created.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    mark_freq : list of float64, optional
         (Default value = '')
         A list of frequencies to mark on the plot with dashed lines.
    z : float64, optional
         (Default value = 0.)
         Redshift by which the frequencies given in 'mark_freq' and 'mark_CO' parameters should be shifted.
         Currently only one redshift can be given for all targets.
    mark_CO : bool, optional
         (Default value = False)
         Mark CO, 13CO, and C18O frequencies on the plot with dashed lines.

    Notes
    -----
    Reports will be grouped by ALMA target names, therefore the same source with many different ALMA names will be
    treated as individual unique targets (e.g. TW_Hya, TW Hya, twhya).

    """
    source_list = observations['target_name'].unique()
    for s, source in enumerate(source_list):
        detections = observations[observations['target_name'] == source]
        plot_overview(detections, mark_freq=mark_freq, z=z, mark_CO=mark_CO, showfig=False, savefig=source)


##############################################
# Downloading Data
##############################################

def download_data(observations, fitsonly=False, dryrun=False, print_urls=False, filename_must_include='',
                  location='./data'):
    """
    Download ALMA data from the archive to a location on the local machine.

    Parameters
    ----------
    observations : pandas.DataFrame
         This is likely the output of e.g. 'conesearch', 'target', 'catalog', & 'keysearch' functions.
    fitsonly : bool, optional
         (Default value = False)
         Download individual fits files only (fitsonly=True). This option will not download the raw data
         (e.g. 'asdm' files), weblogs, or README files.
    dryrun : bool, optional
         (Default value = False)
         Allow the user to do a test run to check the size and number of files to download without actually
         downloading the data (dryrun=True). To download the data, set dryrun=False.
    print_urls : bool, optional
         (Default value = False)
         Write the list of urls to be downloaded from the archive to the terminal.
    filename_must_include : list of str, optional
         (Default value = '')
         A list of strings the user wants to be contained in the url filename. This is useful to restrict the
         download further, for example, to data that have been primary beam corrected ('.pbcor') or that have
         the science target or calibrators (by including their names). The choice is largely dependent on the
         cycle and type of reduction that was performed and data products that exist on the archive as a result.
         In most recent cycles, the science target can be filtered out with the flag '_sci' or its ALMA target name.
    location : str, optional
         (Default value = ./data)
         directory where the downloaded data should be placed.

    """
    print("================================")
    # we use astroquery to download data
    myAlma = Alma()
    default_location = './data'
    myAlma.cache_location = default_location
    # catch the case where the DataFrame is empty.
    try:
        if any(observations['data_rights'] == 'Proprietary'):
            print("Warning: some of the data you are trying to download are still in the proprietary period and are "
                  "not publicly available yet.")
            observations = observations[observations['data_rights'] == 'Public']
        uids_list = observations['member_ous_uid'].unique()
        # when len(uids_list) == 0, it's because the DataFrame included only proprietary data and we removed them in
        # the above if statement, so the DataFrame is now empty
        if len(uids_list) == 0:
            print("No data to download. Check the input DataFrame. It is likely that your query results include only "
                  "proprietary data which cannot be freely downloaded.")
            return
    # this is the case where the query had no results to begin with.
    except TypeError:
        print("No data to download. Check the input DataFrame.")
        return
    # change download location if specified by user, else the location will be the astrquery cache location
    if location != default_location:
        if os.path.isdir(location):
            myAlma.cache_location = location
        else:
            print("{} is not a directory. The download location will be set to {}".format(location, default_location))
            myAlma.cache_location = default_location
    if fitsonly:
        data_table = Alma.get_data_info(uids_list, expand_tarfiles=True)
        # filter the data_table and keep only rows with "fits" in 'access_url' and the strings provided by user
        # in 'filename_must_include' parameter
        dl_table = data_table[[i for i, v in enumerate(data_table['access_url']) if v.endswith(".fits") and
                               all(i in v for i in filename_must_include)]]
        dl_link_list = dl_table['access_url'].tolist()
        # keep track of the download size and number of files to download
        dl_size = dl_table['content_length'].sum() / 1E9
        dl_files = len(dl_table)
        if dryrun:
            print("This is a dryrun. To begin download, set dryrun=False.")
            print("================================")
        else:
            print("Starting download. Please wait...")
            print("================================")
            myAlma.download_files(dl_link_list, cache=True)
    else:
        data_table = Alma.get_data_info(uids_list, expand_tarfiles=False)
        dl_link_list = data_table['access_url'].tolist()
        # keep track of the download size and number of files to download
        dl_size = data_table['content_length'].sum() / 1E9
        dl_files = len(data_table)
        if dryrun:
            print("This is a dryrun. To begin download, set dryrun=False.")
            print("================================")
        else:
            print("Starting download. Please wait...")
            print("================================")
            myAlma.retrieve_data_from_uid(uids_list, cache=True)
    print("Download location = {}".format(myAlma.cache_location))
    print("Total number of Member OUSs to download = {}".format(len(uids_list)))
    print("Selected Member OUSs: {}".format(uids_list.tolist()))
    print("Number of files to download = {}".format(dl_files))
    print("Needed disk space = {:.1f} GB".format(dl_size))
    if print_urls:
        print("File URLs to download = {}".format("\n".join(dl_link_list)))
    print("--------------------------------")

####################################################################
