"""
Script to create full data release for PHANGS-HST cluster catalogs
"""

from data_release_routines import DataReleaseRoutines

try:
    import data_release_config
except ImportError:
    raise ImportError('No data_access_config.py file found. This file needs to be created to specify data paths to '
                      'internal release etc. An example is shown in the file data_access_config_example.py ')

_author = 'Daniel Maschmann'
_contact = 'dmaschmann(at)stsci(dot)edu'
_data_release = 5
_catalog_release = 2

# get local data path structure
config_dict = data_release_config.config_dict

# Create object of DataReleaseRoutine class and give it all variables from the config_dict
data_release_ojt = DataReleaseRoutines(**config_dict)

# create the final data release (all the catalogues)
data_release_ojt.create_final_data_release()
