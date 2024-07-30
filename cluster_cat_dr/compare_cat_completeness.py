"""
Script to create full data release for PHANGS-HST cluster catalogs
"""

from astropy.io import fits
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
data_release_ojt_new = DataReleaseRoutines(**config_dict)
data_release_ojt_old = DataReleaseRoutines(**config_dict)

data_release_ojt_old.cat_ver = 'v2'

target_list = ['ic1954', 'ic5332', 'ngc0628e', 'ngc0628c', 'ngc0685', 'ngc1087', 'ngc1097', 'ngc1300',
                            'ngc1317', 'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1559', 'ngc1566', 'ngc1672',
                            'ngc1792', 'ngc2775', 'ngc2835', 'ngc2903', 'ngc3351', 'ngc3621', 'ngc3627', 'ngc4254',
                            'ngc4298', 'ngc4303', 'ngc4321', 'ngc4535', 'ngc4536', 'ngc4548', 'ngc4569', 'ngc4571',
                            'ngc4654', 'ngc4689', 'ngc4826', 'ngc5068', 'ngc5248', 'ngc6744', 'ngc7496']





for target in target_list:

    # get daves catalog
    daves_selection = data_release_ojt_new.get_ir_cat(target, table_number=1, classify='ml',
                                                         cl_class='class12')

    cat_name_new = data_release_ojt_new.get_data_release_table_name(target=target, classify='ml',
                                                            cl_class='class12', table_type='obs',cat_ver='v1', sed_type=None)

    cat_name_old = data_release_ojt_old.get_data_release_table_name(target=target, classify='ml',
                                                            cl_class='class12', table_type='obs',cat_ver='v2', sed_type=None)

    data_new_cat = fits.open('/home/benutzer/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2_IR5/catalogs/' + cat_name_new)[1].data
    # get old catalog
    data_old_cat = fits.open('/home/benutzer/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2/catalogs/' + cat_name_old)[1].data


    print(target, len(daves_selection), len(data_new_cat),  len(data_old_cat),  len(data_new_cat) - len(data_old_cat))


