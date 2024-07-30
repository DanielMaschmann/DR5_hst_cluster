"""
specification configurations and data structures for the local machine
"""

""" all information and configurations which is specific for the local machine must be gathered in this dictionary """
config_dict = {

    # name of produced DR
    'hst_cc_ver': 'IR5',

    # names of catalogs which will be included
    'sed_cat_1': 'opt_and_info_ground_based_ha',
    'sed_cat_2': 'opt_and_info_hst_ha',
    'sed_cat_3': 'opt_and_hst_ha',
    'sed_cat_4': 'opt_and_hst_ha_and_nir',

    # path to raw data tables
    'path2sed_cat_1': '/home/benutzer/data/PHANGS_products/HST_catalogs/'
                      'SEDfix_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_inclusiveGCclass_Jun21',
    'path2sed_cat_2': '/home/benutzer/data/PHANGS_products/HST_catalogs/'
                      'SEDfix_NewModelsHSTHaUnionHaFLAG11pc_inclusiveGCcc_inclusiveGCclass_Jun21',
    'path2sed_cat_3': '/media/benutzer/Extreme Pro/data/phangs_data_products/kiana_sed_fit/NUV_optical_NIR_final_fits/',
    'path2sed_cat_4': '/media/benutzer/Extreme Pro/data/phangs_data_products/kiana_sed_fit/NUV_optical_NIR_final_fits/',

    # set path to column tables
    'path2col_name_tab_1': '/home/benutzer/Documents/projects/DR5_hst_cluster/cluster_cat_dr/column_data/'
                           'DR5_sed_column_name_table - col_names_dr4_cr2.csv',
    'path2col_name_tab_2': '/home/benutzer/Documents/projects/DR5_hst_cluster/cluster_cat_dr/column_data/'
                           'DR5_sed_column_name_table - col_names_dr4_cr2.csv',
    'path2col_name_tab_3': '',
    'path2col_name_tab_4': '',

    # which galaxies will be included in the
    'target_list_table_1': ['ic1954', 'ic5332', 'ngc0628e', 'ngc0628c', 'ngc0685', 'ngc1087', 'ngc1097', 'ngc1300',
                            'ngc1317', 'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1559', 'ngc1566', 'ngc1672',
                            'ngc1792', 'ngc2775', 'ngc2835', 'ngc2903', 'ngc3351', 'ngc3621', 'ngc3627', 'ngc4254',
                            'ngc4298', 'ngc4303', 'ngc4321', 'ngc4535', 'ngc4536', 'ngc4548', 'ngc4569', 'ngc4571',
                            'ngc4654', 'ngc4689', 'ngc4826', 'ngc5068', 'ngc5248', 'ngc6744', 'ngc7496'],
    'target_list_table_2': ['ic5332', 'ngc0628e', 'ngc0628c', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433',
                            'ngc1512', 'ngc1566', 'ngc1672', 'ngc3351', 'ngc3627', 'ngc4254', 'ngc4303', 'ngc4321',
                            'ngc5068', 'ngc7496'],
    'target_list_table_3': ['ic5332', 'ngc0628e', 'ngc0628c', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433',
                            'ngc1512', 'ngc1566', 'ngc1672', 'ngc3351', 'ngc3627', 'ngc4254', 'ngc4303', 'ngc4321',
                            'ngc5068', 'ngc7496'],
    'target_list_table_4': ['ic5332', 'ngc0628e', 'ngc0628c', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433',
                            'ngc1512', 'ngc1566', 'ngc1672', 'ngc3351', 'ngc3627', 'ngc4254', 'ngc4303', 'ngc4321',
                            'ngc5068', 'ngc7496'],

    'create_candidate_table_1': True,
    'create_candidate_table_2': True,
    'create_candidate_table_3': False,
    'create_candidate_table_4': False,

    # parameters to create data out put
    # note that for DR 4 CR2 we will create Tables 1 and 2 and for DR5 CR1 tables 3 and 4.
    'list_include_tables': [1, 2],
    # Data output
    # data release
    'data_release': '4',
    # catalog release
    'catalog_release': '2',
    # catalog version
    'cat_ver': 'v1',
    # path to Data Release 4, Catalog Release 3 at MAST
    'catalog_output_path': '/home/benutzer/data/PHANGS_products/HST_catalogs',

    # Identified artifact in the internal release
    # path to tables/catalogs of artifacts
    'path2artifact': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/AR1',
    # a re-evaluation of questionable artefacts.
    'path2questionable_artifacts': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/'
                                   'Questionable Artifacts - questionable_artifacts.csv',
    # masks for diffraction spikes
    'path2diffraction_spike_masks': '/home/benutzer/data/PHANGS_products/HST_catalogs/Artifact_Removal/'
                                    'diffraction_spikes/fits',
    # limit for V-I color and concentration index
    'v_i_color_lim': 2.0,
    'ci_lim': 1.45,

    # flag if artifact removal step should be done
    'artifact_removal_flag': True,
    # flag to raise error if a file is not found. This can be handy for early stages in the development as not every
    # target has yet an artifact catalog. For the final version, however, this should be set to True
    'artifact_rais_file_not_found_flag': True,

    # additional
    'existing_artifact_removal_flag': True,

}
