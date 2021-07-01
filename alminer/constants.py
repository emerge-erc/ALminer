"""Global variables & constants to be used in alminer.py"""

# Lists and dictionaries for analysis and plotting purposes
band_names = ["Band 3", "Band 4", "Band 5", "Band 6", "Band 7", "Band 8, 9, 10"]

band_color = {"Band 3": "#BDD9BF", "Band 4": "#929084", "Band 5": "#FFC857", "Band 6": "#A997DF",
              "Band 7": "#E5323B", "Band 8, 9, 10": "#2E4052"}

band_min_freq = {'Band 3': 84., 'Band 4': 125., 'Band 5': 163., 'Band 6': 211., 'Band 7': 275., 'Band 8, 9, 10': 373.}

band_max_freq = {'Band 3': 116., 'Band 4': 163., 'Band 5': 211., 'Band 6': 275., 'Band 7': 373., 'Band 8, 9, 10': 950.}

# CO, 13CO, C180 lines covered in all ALMA bands
CO_line_names = ["CO (1-0)", "CO (2-1)", "CO (3-2)", "CO (4-3)", "CO (5-4)", "CO (6-5)",
                 "CO (7-6)", "CO (8-7)", "13CO (1-0)", "13CO (2-1)", "13CO (3-2)",
                 "13CO (4-3)",  "13CO (5-4)", "13CO (6-5)", "13CO (7-6)", "13CO (8-7)",
                 "C18O (1-0)", "C18O (2-1)", "C18O (3-2)", "C18O (4-3)", "C18O (5-4)", "C18O (6-5)",
                 "C18O (7-6)", "C18O (8-7)"]

CO_line_freq = {"CO (1-0)": 115.27120180, "CO (2-1)": 230.53800000, "CO (3-2)": 345.79598990,
                "CO (4-3)": 461.0407682, "CO (5-4)": 576.2679305, "CO (6-5)": 691.4730763,
                "CO (7-6)": 806.6518060, "CO (8-7)": 921.7997000,
                "13CO (1-0)": 110.20132180, "13CO (2-1)": 220.39861950, "13CO (3-2)": 330.58786710,
                "13CO (4-3)": 440.7651735, "13CO (5-4)": 550.9262851, "13CO (6-5)": 661.0672766,
                "13CO (7-6)": 771.1841250, "13CO (8-7)": 881.2728080,
                "C18O (1-0)": 109.78217340, "C18O (2-1)": 219.56035410, "C18O (3-2)": 329.33055250,
                "C18O (4-3)": 439.0887658, "C18O (5-4)": 548.8310055, "C18O (6-5)": 658.5532782,
                "C18O (7-6)": 768.2515933, "C18O (8-7)": 877.9219553}

CO_line_ha = {"CO (1-0)": "left", "CO (2-1)": "left", "CO (3-2)": "left", "CO (4-3)": "left", "CO (5-4)": "left",
              "CO (6-5)": "left", "CO (7-6)": "left", "CO (8-7)": "left", "13CO (1-0)": "left", "13CO (2-1)": "left",
              "13CO (3-2)": "left", "13CO (4-3)": "left",  "13CO (5-4)": "left", "13CO (6-5)": "left",
              "13CO (7-6)": "left",  "13CO (8-7)": "left", "C18O (1-0)": "right", "C18O (2-1)": "right",
              "C18O (3-2)": "right", "C18O (4-3)": "right", "C18O (5-4)": "right", "C18O (6-5)": "right",
              "C18O (7-6)": "right", "C18O (8-7)": "right"}

CO_line_label = {"CO (1-0)": r'$\mathregular{\,CO\,(1-0)\,\,}$', "CO (2-1)": r'$\mathregular{\,CO\,(2-1)\,\,}$',
                 "CO (3-2)": r'$\mathregular{\,CO\,(3-2)\,\,}$', "CO (4-3)": r'$\mathregular{\,CO\,(4-3)\,\,}$',
                 "CO (5-4)": r'$\mathregular{\,CO\,(5-4)\,\,}$', "CO (6-5)": r'$\mathregular{\,CO\,(6-5)\,\,}$',
                 "CO (7-6)": r'$\mathregular{\,CO\,(7-6)\,\,}$', "CO (8-7)": r'$\mathregular{\,CO\,(8-7)\,\,}$',
                 "13CO (1-0)": r'$\mathregular{\,^{13}CO (1-0)}$',
                 "13CO (2-1)": r'$\mathregular{\,^{13}CO\,(2-1)\,\,}$',
                 "13CO (3-2)": r'$\mathregular{\,^{13}CO\,(3-2)\,\,}$',
                 "13CO (4-3)": r'$\mathregular{\,^{13}CO\,(4-3)\,\,}$',
                 "13CO (5-4)": r'$\mathregular{\,^{13}CO\,(5-4)\,\,}$',
                 "13CO (6-5)": r'$\mathregular{\,^{13}CO\,(6-5)\,\,}$',
                 "13CO (7-6)": r'$\mathregular{\,^{13}CO\,(7-6)\,\,}$',
                 "13CO (8-7)": r'$\mathregular{\,^{13}CO\,(8-7)\,\,}$',
                 "C18O (1-0)": r'$\mathregular{\,C^{18}O\,(1-0)\,\,}$',
                 "C18O (2-1)": r'$\mathregular{\,C^{18}O\,(2-1)\,\,}$',
                 "C18O (3-2)": r'$\mathregular{\,C^{18}O\,(3-2)\,\,}$',
                 "C18O (4-3)": r'$\mathregular{\,C^{18}O\,(4-3)\,\,}$',
                 "C18O (5-4)": r'$\mathregular{\,C^{18}O\,(5-4)\,\,}$',
                 "C18O (6-5)": r'$\mathregular{\,C^{18}O\,(6-5)\,\,}$',
                 "C18O (7-6)": r'$\mathregular{\,C^{18}O\,(7-6)\,\,}$',
                 "C18O (8-7)": r'$\mathregular{\,C^{18}O\,(8-7)\,\,}$'}

# Define all possible keywords from TAP and their types
VALID_KEYWORDS_STR = ('obs_publisher_did', 'obs_collection', 'facility_name', 'instrument_name', 'obs_id',
                      'dataproduct_type', 'target_name', 's_region', 'pol_states', 'o_ucd', 'band_list',
                      'authors', 'pub_abstract', 'proposal_abstract', 'schedblock_name', 'proposal_authors',
                      'group_ous_uid', 'member_ous_uid', 'asdm_uid', 'obs_title', 'type', 'scan_intent',
                      'science_observation', 'antenna_arrays', 'is_mosaic', 'obs_release_date', 'frequency_support',
                      'obs_creator_name', 'pub_title', 'first_author', 'qa2_passed', 'bib_reference',
                      'science_keyword', 'scientific_category', 'lastModified', 'access_url', 'access_format',
                      'proposal_id', 'data_rights')
VALID_KEYWORDS_DOU = ('gal_longitude', 'gal_latitude', 's_ra', 's_dec', 's_fov', 's_resolution', 't_min', 't_max',
                      't_exptime', 't_resolution', 'em_min', 'em_max', 'em_res_power', 'em_resolution',
                      'sensitivity_10kms', 'cont_sensitivity_bandwidth', 'spatial_scale_max', 'bandwidth',
                      'spatial_resolution', 'frequency', 'velocity_resolution', 'pwv')
VALID_KEYWORDS_INT = ('calib_level', 'publication_year')

# Define new columns we manually add to the dataframe
NEW_COLUMNS = ['Obs', 'project_code', 'ALMA_source_name', 'RAJ2000', 'DEJ2000', 'ang_res_arcsec', 'min_freq_GHz',
               'max_freq_GHz', 'central_freq_GHz', 'bandwidth_GHz', 'freq_res_kHz', 'vel_res_kms', 'LAS_arcsec',
               'FoV_arcsec', 'cont_sens_bandwidth', 'line_sens_10kms', 'line_sens_native', 'MOUS_id']

# setting the column types for the manually added columns
COLUMN_TYPES = {'Obs': 'Int64', 'project_code': str, 'ALMA_source_name': str, 'RAJ2000': 'float64',
                'DEJ2000': 'float64', 'ang_res_arcsec': 'float64', 'min_freq_GHz': 'float64',
                'max_freq_GHz': 'float64',
                'central_freq_GHz': 'float64', 'bandwidth_GHz': 'float64', 'freq_res_kHz': 'float64',
                'vel_res_kms': 'float64', 'LAS_arcsec': 'float64', 'FoV_arcsec': 'float64',
                'cont_sens_bandwidth': 'float64', 'line_sens_10kms': 'float64', 'line_sens_native': 'float64',
                'MOUS_id': str}

# setting the column types for the TAP columns
COLUMN_TYPES.update({k: 'float64' for k in VALID_KEYWORDS_DOU})
COLUMN_TYPES.update({k: 'Int64' for k in VALID_KEYWORDS_INT})
COLUMN_TYPES.update({k: str for k in VALID_KEYWORDS_STR})
COLUMN_TYPES.update({'publication_year': object})
