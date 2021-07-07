"""
ALminer: ALMA archive mining and visualization toolkit
=======
A package for mining the Atacama Large Millimeter/submillimeter
Array (ALMA) data archive and visualizing the queried observations.
"""
from .alminer import catalog, CO_lines, conesearch, download_data, explore, filter_results, get_description, \
    get_units, get_info, keysearch, line_coverage, plot_line_overview, plot_overview, plot_observations, plot_bands, \
    plot_sky, run_query, summary, target, save_source_reports, save_table
__all__ = ["catalog", "CO_lines", "conesearch", "download_data", "explore", "filter_results",
           "get_description", "get_units", "get_info", "keysearch", "line_coverage", "plot_line_overview",
           "plot_overview", "plot_observations", "plot_bands", "plot_sky", "run_query", "summary",
           "target", "save_source_reports", "save_table"]
try:
    from .version import __version__
except ImportError:
    __version__ = " "
