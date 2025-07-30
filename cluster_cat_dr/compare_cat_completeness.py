"""
Script to create full data release for PHANGS-HST cluster catalogs
"""

from astropy.io import fits
from data_release_routines import DataReleaseRoutines
from phangs_data_access import helper_func
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



    cat_name_hum12 = data_release_ojt_new.get_data_release_table_name(target=target, classify='human', cl_class='class12', table_type='obs',cat_ver='v2', sed_type=None)
    cat_name_hum3 = data_release_ojt_new.get_data_release_table_name(target=target, classify='human', cl_class='class3', table_type='obs',cat_ver='v2', sed_type=None)
    cat_name_ml12 = data_release_ojt_new.get_data_release_table_name(target=target, classify='ml', cl_class='class12', table_type='obs',cat_ver='v2', sed_type=None)
    cat_name_ml3 = data_release_ojt_new.get_data_release_table_name(target=target, classify='ml', cl_class='class3', table_type='obs',cat_ver='v2', sed_type=None)

    data_old_hum12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2/catalogs/' + cat_name_hum12)[1].data
    data_old_hum3 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2/catalogs/' + cat_name_hum3)[1].data
    data_old_ml12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2/catalogs/' + cat_name_ml12)[1].data
    data_old_ml3 = fits.open( '/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr2/catalogs/' + cat_name_ml3)[1].data

    cat_name_hum12 = data_release_ojt_new.get_data_release_table_name(target=target, classify='human', cl_class='class12', table_type='obs',cat_ver='v1', sed_type=None)
    cat_name_hum3 = data_release_ojt_new.get_data_release_table_name(target=target, classify='human', cl_class='class3', table_type='obs',cat_ver='v1', sed_type=None)
    cat_name_ml12 = data_release_ojt_new.get_data_release_table_name(target=target, classify='ml', cl_class='class12', table_type='obs',cat_ver='v1', sed_type=None)
    cat_name_ml3 = data_release_ojt_new.get_data_release_table_name(target=target, classify='ml', cl_class='class3', table_type='obs',cat_ver='v1', sed_type=None)

    data_new_hum12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr5_cr1_IR5/catalogs/' + cat_name_hum12)[1].data
    data_new_hum3 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr5_cr1_IR5/catalogs/' + cat_name_hum3)[1].data
    data_new_ml12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr5_cr1_IR5/catalogs/' + cat_name_ml12)[1].data
    data_new_ml3 = fits.open( '/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr5_cr1_IR5/catalogs/' + cat_name_ml3)[1].data

    cat_name_no_cut_hum12 = 'SEDfix_PHANGS_IR4_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21_phangs_hst_v1p2_human_class12.fits' % helper_func.FileTools.target_names_no_zeros(target=target)
    cat_name_no_cut_hum3 = 'SEDfix_PHANGS_IR4_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21_phangs_hst_v1p2_human_class3.fits' % helper_func.FileTools.target_names_no_zeros(target=target)
    cat_name_no_cut_ml12 = 'SEDfix_PHANGS_IR4_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21_phangs_hst_v1p2_ml_class12.fits' % helper_func.FileTools.target_names_no_zeros(target=target)
    cat_name_no_cut_ml3 = 'SEDfix_PHANGS_IR4_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21_phangs_hst_v1p2_ml_class3.fits' % helper_func.FileTools.target_names_no_zeros(target=target)

    data_no_cut_hum12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/SEDfix_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21/' + cat_name_no_cut_hum12)[1].data
    data_no_cut_hum3 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/SEDfix_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21/' + cat_name_no_cut_hum3)[1].data
    data_no_cut_ml12 = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/SEDfix_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21/' + cat_name_no_cut_ml12)[1].data
    data_no_cut_ml3 = fits.open( '/Users/dmaschmann/data/PHANGS_products/HST_catalogs/SEDfix_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21/' + cat_name_no_cut_ml3)[1].data


    print(target,
          'HUMC12 ', len(data_old_hum12), len(data_new_hum12), len(data_new_hum12) - len(data_old_hum12),'  |  ',
          'HUMC3 ', len(data_old_hum3), len(data_new_hum3), len(data_new_hum3) - len(data_old_hum3),'  |  ',
          'MLC12 ', len(data_old_ml12), len(data_new_ml12), len(data_new_ml12) - len(data_old_ml12),'  |  ',
          'MLC3 ', len(data_old_ml3), len(data_new_ml3), len(data_new_ml3) - len(data_old_ml3),'  |  ',

          )

    # get old catalog
    # data_old_cat = fits.open('/Users/dmaschmann/data/PHANGS_products/HST_catalogs/phangs_hst_cc_dr4_cr3_ground_based_ha/catalogs/' + cat_name_old)[1].data


    # print(target, len(daves_selection), len(data_new_cat),  len(data_old_cat),  len(data_new_cat) - len(data_old_cat))


