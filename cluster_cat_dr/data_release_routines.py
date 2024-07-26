"""
All production routines for the final PHANGS-HST cluster catalog data release will be gathered here
"""

import os
from datetime import datetime
from pathlib import Path
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from astropy.table import Column, hstack, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import ConvexHull

from phangs_data_access import helper_func



class DataReleaseRoutines:
    def __init__(self,
                 hst_cc_ver,
                 sed_cat_1,
                 sed_cat_2,
                 sed_cat_3,
                 sed_cat_4,
                 path2sed_cat_1,
                 path2sed_cat_2,
                 path2sed_cat_3,
                 path2sed_cat_4,
                 path2col_name_tab_1,
                 path2col_name_tab_2,
                 path2col_name_tab_3,
                 path2col_name_tab_4,
                 target_list_table_1,
                 target_list_table_2,
                 target_list_table_3,
                 target_list_table_4,
                 create_candidate_table_1,
                 create_candidate_table_2,
                 create_candidate_table_3,
                 create_candidate_table_4,
                 list_include_tables,
                 data_release,
                 catalog_release,
                 cat_ver,
                 catalog_output_path,
                 path2artifact,
                 path2questionable_artifacts,
                 path2diffraction_spike_masks,
                 v_i_color_lim,
                 ci_lim,
                 artifact_removal_flag,
                 artifact_rais_file_not_found_flag,
                 existing_artifact_removal_flag,
                 ):

        # add input keys to attributes
        self.hst_cc_ver = hst_cc_ver

        self.sed_cat_1 = sed_cat_1
        self.sed_cat_2 = sed_cat_2
        self.sed_cat_3 = sed_cat_3
        self.sed_cat_4 = sed_cat_4
        self.path2sed_cat_1 = path2sed_cat_1
        self.path2sed_cat_2 = path2sed_cat_2
        self.path2sed_cat_3 = path2sed_cat_3
        self.path2sed_cat_4 = path2sed_cat_4
        self.path2col_name_tab_1 = path2col_name_tab_1
        self.path2col_name_tab_2 = path2col_name_tab_2
        self.path2col_name_tab_3 = path2col_name_tab_3
        self.path2col_name_tab_4 = path2col_name_tab_4
        self.target_list_table_1 = target_list_table_1
        self.target_list_table_2 = target_list_table_2
        self.target_list_table_3 = target_list_table_3
        self.target_list_table_4 = target_list_table_4
        self.create_candidate_table_1 = create_candidate_table_1
        self.create_candidate_table_2 = create_candidate_table_2
        self.create_candidate_table_3 = create_candidate_table_3
        self.create_candidate_table_4 = create_candidate_table_4

        self.list_include_tables = list_include_tables

        self.data_release = data_release
        self.catalog_release = catalog_release
        self.cat_ver = cat_ver

        # add version folder to catalog_output_path
        self.catalog_output_path = catalog_output_path + '/phangs_hst_cc_dr%s_cr%s_%s' % (data_release, catalog_release,
                                                                                          hst_cc_ver)
        # artifacts
        self.path2artifact = path2artifact
        self.path2questionable_artifacts = path2questionable_artifacts
        self.path2diffraction_spike_masks = path2diffraction_spike_masks
        self.v_i_color_lim = v_i_color_lim
        self.ci_lim = ci_lim
        self.artifact_removal_flag = artifact_removal_flag
        self.artifact_rais_file_not_found_flag = artifact_rais_file_not_found_flag
        self.existing_artifact_removal_flag = existing_artifact_removal_flag

        self.questionable_artefact_table = None

        # load constructor of parent class
        super().__init__()

    def create_final_data_release(self):
        """
        Function to run all steps to create the PHANGs-HST cluster catalog data release
        """
        # create catalogs
        # loop over all targets for which a cluster catalog exists

        # loop over targets which have to be created
        for table_index in self.list_include_tables:

            print('table_index ', table_index, ' Creating Table ', getattr(self, 'sed_cat_%i' % table_index))
            # loop over target list
            target_list = getattr(self, 'target_list_table_%i' % table_index)
            for target in target_list:
                print('target ', target)

                # get artifact removal table created by Chris and Lucious
                table_artifact = self.get_artifact_cat(target=target)
                # get table with second classification by BCW.
                # This table is only treating the objects Chris and Lucious were not sure about
                table_re_classified = self.get_second_classify_cat(target=target)
                # get table for diffraction spikes
                # This is a little bit hacky...
                # First, We just used a quick way to find an exception for the NGC1512 where we only have a central map
                # second, for NGC 628 there are two fields (center and east).
                # The east field is called ngc628
                # The center field is called ngc628c
                # in the central case, the maps would be loaded into `diffraction_spike_mask`
                # in the east case the maps would be loaded into `diffraction_spike_mask_2`
                # of course there are better ways. But this is working and loading the different band names can be
                # tricky when there is a c or e in the name
                if target == 'ngc1512':
                    diffraction_spike_mask, diffraction_spike_wcs = self.load_diffraction_spike_masks(target=target, target_str='ngc1512c')
                else:
                    diffraction_spike_mask, diffraction_spike_wcs = self.load_diffraction_spike_masks(target=target, target_str=target)
                if (target == 'ngc0628e') | (target == 'ngc0628c'):
                    diffraction_spike_mask_2, diffraction_spike_wcs_2 = (
                        self.load_diffraction_spike_masks(target='ngc0628', target_str='ngc0628'))
                else:
                    diffraction_spike_mask_2 = None
                    diffraction_spike_wcs_2 = None

                # for table 1 and 2, we have to first create the candidate tables
                # this will be the basis for the obs and sed catalogs.
                #
                # for tables 3 and 4 we first create the human class 1+2 catalog created for table 1 and 2
                # we then create a blank catalog and fill the new data in a row matched fashion in.
                #
                if table_index in [1, 2]:
                    # load the raw candidate table
                    candidate_table_ir = self.get_ir_cat(target, table_number=table_index, classify=None,
                                                         cl_class='candidates')
                    # find cross-matching artifact
                    for artifact_index in range(len(table_artifact)):
                        artifact_in_table_ir = ((candidate_table_ir['ID_PHANGS_CLUSTERS_v1p2'] ==
                                                 table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]) &
                                                (candidate_table_ir['PHANGS_X'] ==
                                                 table_artifact['PHANGS_X'][artifact_index]) &
                                                (candidate_table_ir['PHANGS_Y'] ==
                                                 table_artifact['PHANGS_Y'][artifact_index]) &
                                                (candidate_table_ir['PHANGS_RA'] ==
                                                 table_artifact['PHANGS_RA'][artifact_index]) &
                                                (candidate_table_ir['PHANGS_DEC'] ==
                                                 table_artifact['PHANGS_DEC'][artifact_index]))
                        artifact_in_table_re_classified = \
                            ((table_re_classified['ID_PHANGS_CLUSTERS_v1p2'] ==
                              table_artifact['ID_PHANGS_CLUSTERS_v1p2'][artifact_index]))
                        # Update human classification:
                        #
                        # If this row was already classified by human (BCW) we do nothing
                        if np.invert(np.isnan(candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN']
                                              [artifact_in_table_ir][0])):
                            continue
                        # For the classification there have been two steps:
                        # step 1) Chris and Lucious classified artifacts.
                        #         If They were unable to classify the artifact they put in -999
                        # step 2) All the artifacts not classified by Chris and Lucious (-999) were then re-classified
                        #         by BCW.
                        #
                        # Change the human classification according to Chris and Lucious
                        if table_artifact['NEW_CLASS'][artifact_index] != -999:
                            candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir] = (
                                table_artifact)['NEW_CLASS'][artifact_index]
                        # Change the re-classified according to BCW
                        else:
                            # If there is no data in the column, the entry will be masked. In this case we move on.
                            if not np.ma.isMaskedArray(table_re_classified['BCW_estimate']
                                                       [artifact_in_table_re_classified][0]):
                                # In some cases a -999 is also indicating no classification.
                                if int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0]) != -999:
                                    # Updating the human classification.
                                    candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][artifact_in_table_ir] = \
                                        int(table_re_classified['BCW_estimate'][artifact_in_table_re_classified][0])
                    # Some clusters are situated in diffraction spikes of bright foreground stars.
                    # These objects have bad photometry and are therefore flagged as an artifact (artifact code 8)
                    if diffraction_spike_mask is not None:
                        # loading the RA and DEC position of the clusters
                        ra_pixel_coords = candidate_table_ir['PHANGS_RA']
                        dec_pixel_coords = candidate_table_ir['PHANGS_DEC']
                        pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                        # transforming the position into pixel coordinates using the WCS of the diffraction maps
                        pos_pix = diffraction_spike_wcs.world_to_pixel(pos)
                        # computing the integer pixel index of the cluster position
                        x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                        y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                        # Creating a mask of clusters which are situated in the diffraction spike maps
                        # They need to be inside the diffraction spike image (pixel_coords > 0) and
                        # pixel_coords < diffraction_spike_mask.shape[0]
                        mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                                    (x_pixel_coords < diffraction_spike_mask.shape[0]) &
                                                    (y_pixel_coords < diffraction_spike_mask.shape[1]))
                        # Masking the diffraction spike.
                        # The map is a boolean map so 0 indicated no diffraction spike and 1 means diffraction spike
                        artifact_in_diffraction_spike = (
                                diffraction_spike_mask[y_pixel_coords[mask_covered_coordinates],
                                                       x_pixel_coords[mask_covered_coordinates]] > 0)
                        # Updating the human classification with the diffraction spike artifact code
                        candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][mask_covered_coordinates][artifact_in_diffraction_spike] = 8
                    # In some cases there are two diffraction spike maps.
                    # This is due to the fact that two different bands were evaluated.
                    # we now repeat the exact same steps again for the second map if available.
                    if diffraction_spike_mask_2 is not None:
                        ra_pixel_coords = candidate_table_ir['PHANGS_RA']
                        dec_pixel_coords = candidate_table_ir['PHANGS_DEC']
                        pos = SkyCoord(ra=ra_pixel_coords, dec=dec_pixel_coords, unit=(u.degree, u.degree), frame='fk5')
                        pos_pix = diffraction_spike_wcs_2.world_to_pixel(pos)
                        x_pixel_coords = np.array(np.rint(pos_pix[0]), dtype=int)
                        y_pixel_coords = np.array(np.rint(pos_pix[1]), dtype=int)
                        mask_covered_coordinates = ((x_pixel_coords > 0) & (y_pixel_coords > 0) &
                                                    (x_pixel_coords < diffraction_spike_mask_2.shape[0]) &
                                                    (y_pixel_coords < diffraction_spike_mask_2.shape[1]))
                        artifact_in_diffraction_spike = (
                                diffraction_spike_mask_2[y_pixel_coords[mask_covered_coordinates],
                                                         x_pixel_coords[mask_covered_coordinates]] > 0)
                        candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][mask_covered_coordinates][artifact_in_diffraction_spike] = 8

                    # update very red stars
                    # These stars are identified by the concentration index and the V-I color.
                    # Their artifact code is 19
                    vi_color = candidate_table_ir['PHANGS_F555W_vega_tot'] - candidate_table_ir['PHANGS_F814W_vega_tot']
                    ci = candidate_table_ir['PHANGS_CI']
                    # the thresholds for the selection are specified in the DR configuration file
                    very_red_star_mask_candidate_table_ir = (vi_color > self.v_i_color_lim) & (ci < self.ci_lim)
                    candidate_table_ir['PHANGS_CLUSTER_CLASS_HUMAN'][very_red_star_mask_candidate_table_ir] = 19
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # !!!!! Save candidate table !!!!!
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # For this step we will first create the candidate table by looping over all columns
                    # In this loop we load the cvs files providing the final column names, data types, units and
                    # descriptions
                    candidate_table = self.create_cand_table(target=target, table=candidate_table_ir)
                    # sort them by increasing Y pixel
                    increase_y_sort = np.argsort(candidate_table['PHANGS_Y'])
                    candidate_table = candidate_table[increase_y_sort]
                    # Create the PHANGS cluster ID column
                    # get ids from all objects which have been classified as a cluster
                    indexes_with_clusters = np.where(
                        # can be human cluster
                        (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 1) |
                        (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 2) |
                        (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 3) |
                        # or ML cluster
                        (((candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 1) |
                          (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 2) |
                          (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 3)) &
                         # but must be no human classified artefact
                         # otherwise this object is no cluster and will be sorted out.
                         # This last argument is not redundant it sorts those ML cluster out which are humanly
                         # classified artifacts.
                         ((candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] < 4) |
                          np.isnan(candidate_table['PHANGS_CLUSTER_CLASS_HUMAN']))))
                    # create empty data for column
                    cluster_id_column_data = np.ones(len(candidate_table), dtype=int) * -999
                    # only give increasing values for clusters
                    cluster_id_column_data[indexes_with_clusters] = np.arange(1, len(indexes_with_clusters[0])+1)
                    # create cluster id column
                    cluster_id_column = Column(data=cluster_id_column_data, name='ID_PHANGS_CLUSTER',
                                               dtype=int,
                                               description='PHANGS cluster ID for each individual object classified as '
                                                           'class 1,2 or 3, ordered by increasing Y pixel coordinate')
                    # create also a simple running index column
                    index_column = Column(data=np.arange(1, len(candidate_table) + 1),
                                          name='INDEX',
                                          dtype=int,
                                          description='Running index from 1 to N for each individual target')
                    # add these columns to the candidate table
                    candidate_table = hstack([index_column, cluster_id_column, candidate_table])

                    # Create an ML class column but correct them for human classified artifacts
                    ml_class_corr_data = deepcopy(candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'])
                    mask_identified_artifacts_cand = candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] > 3
                    ml_class_corr_data[mask_identified_artifacts_cand] = (
                        candidate_table)[mask_identified_artifacts_cand]['PHANGS_CLUSTER_CLASS_HUMAN']
                    ml_class_corr_column = Column(data=ml_class_corr_data,
                                                  name='PHANGS_CLUSTER_CLASS_ML_VGG_CORR',
                                                  dtype=int,
                                                  description='Same classification as in column `cluster_class_ml\' '
                                                              'but corrected for human identified artefacts. '
                                                              'Thus, all object classes identified as artefacts in '
                                                              'column `class_human\' are replaced.')
                    # add the new column direct behind the VGG classification columns
                    index_last_classification = np.where(np.array(candidate_table.colnames) ==
                                                         'PHANGS_CLUSTER_CLASS_ML_VGG_QUAL')
                    candidate_table.add_column(ml_class_corr_column, index=index_last_classification[0][0]+1)

                    # We now blank SED results for objects with bad observational coverage
                    #
                    # This evaluation is based on the following estimators:
                    # a) PHANGS_HALPHA_EVAL_FLAG : a flag that tells us if the cluster could be evaluated with H-alpha imaging
                    # b) HaFLAG : a flag that tells us if the cluster passes the H-alpha SED-TreeFit decision point.
                    #            This in particular means that the cluster is assigned to the YNO or the GENERAL branch.
                    #            These objects have the possibility to be more dust embedded and sometimes lack the
                    #            detection of one or two HST broad_band filters.
                    # c) n_good_HST_bands : int-value counting the HST broad bands with reliable measurements and
                    #                           can be maximally 5.
                    #                           This number is estimated through the band uncertainty:
                    #                           (PHANGS_<BAND>_mJy_ERR < 0) OR (PHANGS_<BAND>_mJy_ERR > 1.e19)
                    #                           Here, negative uncertainty values indicate a non-detection and crazy
                    #                           high values indicate (1e19) indicate bad coverage.
                    #
                    # We flag the following objects:
                    # 1) Objects with PHANGS_HALPHA_EVAL_FLAG == False
                    #    These objects will drop through the SED decision tree and cannot be properly assigned to a grid
                    # 2) Objects with (PHANGS_HALPHA_EVAL_FLAG == True) AND (HaFLAG = True) AND (n_good_HST_bands < 3)
                    #    These objects have evaluated and detected H-alpha measurements meaning they are assigned to the
                    #    YNO grid but only 0, 1 or 2 HST bands which is not enough for a reasonable SED fit evaluation
                    # 3) Objects with (PHANGS_HALPHA_EVAL_FLAG == True) AND (HaFLAG = False) AND (n_good_HST_bands < 4)
                    #    These objects have H-alpha observational coverage but the cluster is not associated with
                    #    H-alpha emission. These clusters are assigned to the GENERAL so, we therefore do not expect
                    #    high dust extinction hence objects with 0, 1, 2, or 3 HST bands are not reliable fittable.

                    # computing number of good HST bands:
                    hst_band_list = helper_func.BandTools.get_hst_obs_band_list(target=target,
                                                                                phangs_hst_version='phangs_hst_1')
                    # create int array with initial 5 good bands
                    n_good_hst_bands = np.ones(len(candidate_table), dtype=int) * 5
                    # loop over HST bands
                    for hst_band in hst_band_list:
                        bad_hst_band = ((candidate_table['PHANGS_%s_mJy_ERR' % hst_band] < 0.) |
                                        (candidate_table['PHANGS_%s_mJy_ERR' % hst_band] > 1.e19))
                        # subtract each bad HST band
                        n_good_hst_bands -= np.array(bad_hst_band, dtype=int)
                    # create mask to drop SED fit results
                    mask_drop_sed_values = (np.invert(candidate_table['PHANGS_HALPHA_EVAL_FLAG']) |
                                            (candidate_table['PHANGS_HALPHA_EVAL_FLAG'] & candidate_table['HaFLAG'] &
                                             (n_good_hst_bands < 3)) |
                                            (candidate_table['PHANGS_HALPHA_EVAL_FLAG'] &
                                             np.invert(candidate_table['HaFLAG']) &
                                             (n_good_hst_bands < 4)))
                    # now we blank the SED fit results for these objects
                    for col_name in candidate_table.colnames:
                        if 'SED' in col_name:
                            candidate_table[col_name][mask_drop_sed_values] = np.nan

                    # now finally add two columns for CCD classification
                    self.add_ccd_region_col(target, candidate_table)

                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # !!!!!!!!! save Candidate table !!!!!!!!!
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    #
                    self.save_catalog(target=target, table=candidate_table, classify=None, cl_class=None,
                                      table_type='cand', sed_type=getattr(self, 'sed_cat_%i' % table_index))

                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # !!!!!!! OBS and SED catalogs !!!!!!!
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    # getting column names for obs and SED tables
                    col_table = ascii.read(getattr(self, 'path2col_name_tab_%i' % table_index), format='csv')
                    final_colum_names = col_table['final_col_name']
                    obs_column_names = final_colum_names[np.array(col_table['obs_tab'], dtype=bool)]
                    sed_column_names = final_colum_names[np.array(col_table['sed_tab'], dtype=bool)]

                    # get the band name of the u-band
                    band_list = helper_func.BandTools.get_hst_obs_band_list(target=target, phangs_hst_version='phangs_hst_1')
                    if 'F435W' in band_list:
                        not_u_band = 'F438W'
                    elif 'F438W' in band_list:
                        not_u_band = 'F435W'
                    else:
                        raise RuntimeError('Something went wrong. For the galaxy ', target, ' There is no U-band available! It must be F435W or F438W ')

                    # now construct obs table
                    obs_table = None
                    for obs_col_name in obs_column_names:
                        if not_u_band in obs_col_name:
                            continue
                        if obs_table is None:
                            obs_table = candidate_table[obs_col_name]
                        else:
                            obs_table = hstack([obs_table, candidate_table[obs_col_name]])

                    # now construct sed table
                    sed_table = None
                    for sed_col_name in sed_column_names:
                        if not_u_band in sed_col_name:
                            continue
                        if sed_table is None:
                            sed_table = candidate_table[sed_col_name]
                        else:
                            sed_table = hstack([sed_table, candidate_table[sed_col_name]])

                    # now save final catalogs
                    # select human and machine learning catalogs and class 1+2 and class 3
                    # human classified catalogs
                    mask_hum_class12 = ((candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 1) |
                                        (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 2))
                    mask_hum_class3 = candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] == 3
                    # note that we also exclude human classified artefacts
                    mask_ml_class12 = (((candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 1) |
                                       (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 2)) &
                                       (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] > 0) &
                                       (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] < 4))
                    mask_ml_class3 = (candidate_table['PHANGS_CLUSTER_CLASS_ML_VGG'] == 3 &
                                      (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] > 0) &
                                      (candidate_table['PHANGS_CLUSTER_CLASS_HUMAN'] < 4))

                    obs_hum_cass12_cat = obs_table[mask_hum_class12]
                    obs_hum_cass3_cat = obs_table[mask_hum_class3]
                    obs_ml_cass12_cat = obs_table[mask_ml_class12]
                    obs_ml_cass3_cat = obs_table[mask_ml_class3]

                    sed_hum_cass12_cat = sed_table[mask_hum_class12]
                    sed_hum_cass3_cat = sed_table[mask_hum_class3]
                    sed_ml_cass12_cat = sed_table[mask_ml_class12]
                    sed_ml_cass3_cat = sed_table[mask_ml_class3]

                    # re sort the INDEX column again
                    obs_hum_cass12_cat['INDEX'] = np.arange(1, len(obs_hum_cass12_cat) + 1)
                    obs_hum_cass3_cat['INDEX'] = np.arange(1, len(obs_hum_cass3_cat) + 1)
                    obs_ml_cass12_cat['INDEX'] = np.arange(1, len(obs_ml_cass12_cat) + 1)
                    obs_ml_cass3_cat['INDEX'] = np.arange(1, len(obs_ml_cass3_cat) + 1)
                    sed_hum_cass12_cat['INDEX'] = np.arange(1, len(sed_hum_cass12_cat) + 1)
                    sed_hum_cass3_cat['INDEX'] = np.arange(1, len(sed_hum_cass3_cat) + 1)
                    sed_ml_cass12_cat['INDEX'] = np.arange(1, len(sed_ml_cass12_cat) + 1)
                    sed_ml_cass3_cat['INDEX'] = np.arange(1, len(sed_ml_cass3_cat) + 1)

                    # saving the catalogs
                    self.save_catalog(target=target, table=obs_hum_cass12_cat,  classify='human', cl_class='class12',
                                      table_type='obs', sed_type=None)
                    self.save_catalog(target=target, table=obs_hum_cass3_cat,  classify='human', cl_class='class3',
                                      table_type='obs', sed_type=None)
                    self.save_catalog(target=target, table=obs_ml_cass12_cat,  classify='ml', cl_class='class12',
                                      table_type='obs', sed_type=None)
                    self.save_catalog(target=target, table=obs_ml_cass3_cat,  classify='ml', cl_class='class3',
                                      table_type='obs', sed_type=None)

                    self.save_catalog(target=target, table=sed_hum_cass12_cat,  classify='human', cl_class='class12',
                                      table_type='sed', sed_type=getattr(self, 'sed_cat_%i' % table_index))
                    self.save_catalog(target=target, table=sed_hum_cass3_cat,  classify='human', cl_class='class3',
                                      table_type='sed', sed_type=getattr(self, 'sed_cat_%i' % table_index))
                    self.save_catalog(target=target, table=sed_ml_cass12_cat,  classify='ml', cl_class='class12',
                                      table_type='sed', sed_type=getattr(self, 'sed_cat_%i' % table_index))
                    self.save_catalog(target=target, table=sed_ml_cass3_cat,  classify='ml', cl_class='class3',
                                      table_type='sed', sed_type=getattr(self, 'sed_cat_%i' % table_index))

    def save_catalog(self, target, table, classify, cl_class, table_type, sed_type):
        """
        Function to save the catalog tables
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified
        table : ``astropy.table.Table``
        classify : str
            cluster classification either `human` or `ml` or None for candidates
        cl_class : str
            cluster classes must be either `class12` or `class3` or None for candidates
        table_type : str
            to distinguish between obs and sed table
        sed_type : str
        Returns
        -------
        table_ir : ``astropy.io.fits.fitsrec.FITS_rec``
            internal release table
        """

        # check if folder already exists
        if not os.path.isdir(self.catalog_output_path):
            os.mkdir(self.catalog_output_path)
        if not os.path.isdir(self.catalog_output_path + '/catalogs'):
            os.mkdir(self.catalog_output_path + '/catalogs')
        #
        # before we save the table we need to make an exception For the galaxy NGC 1512.
        # This galaxy has a companion NGC 1510 which is present in our observations
        # we simply split the tables into two for this object.
        if target == 'ngc1512':
            # getting the pixel positions
            x_pixel = table['PHANGS_X']
            y_pixel = table['PHANGS_Y']
            # we draw an almost arbitrary line get straight line to divide the points
            x1 = 6107
            y1 = 2580
            x2 = 9408
            y2 = 7010
            slope = (y1-y2)/(x1-x2)
            intersect = y1 - slope * x1
            mask_ngc1510_select = y_pixel < slope * x_pixel + intersect
            # applying the mask
            table_ngc1512 = table[np.invert(mask_ngc1510_select)]
            table_ngc1510 = table[mask_ngc1510_select]
            # rename galaxy name column
            table_ngc1510['PHANGS_GALAXY'] = 'ngc1510'

            # get table name
            cat_name_ngc1512 = self.get_data_release_table_name(target=target, classify=classify, cl_class=cl_class,
                                                                table_type=table_type, sed_type=sed_type)
            # get table name
            cat_name_ngc1510 = self.get_data_release_table_name(target='ngc1510', classify=classify, cl_class=cl_class,
                                                                table_type=table_type, sed_type=sed_type)
            # save table
            table_ngc1512.write(Path(self.catalog_output_path + '/catalogs') /
                                          Path(cat_name_ngc1512),
                                          overwrite=True)
            table_ngc1510.write(Path(self.catalog_output_path + '/catalogs') /
                                          Path(cat_name_ngc1510),
                                          overwrite=True)
        else:
            # get table name
            cat_name = self.get_data_release_table_name(target=target, classify=classify, cl_class=cl_class,
                                                        table_type=table_type, sed_type=sed_type)
            # save table
            table.write(Path(self.catalog_output_path + '/catalogs') / Path(cat_name),  overwrite=True)

    def get_ir_cat(self, target, table_number, classify, cl_class):
        """
        Function to get internal release catalog identifier
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table_number : int
        classify : str
            cluster classification either `human` or `ml`
        cl_class : str
            cluster classes must be either `class12` or `class3`

        Returns
        -------
        table_ir : ``astropy.io.fits.fitsrec.FITS_rec``
            internal release table
        """

        path2table = Path(getattr(self, 'path2sed_cat_%i' % table_number))

        if table_number == 1:
            if cl_class == 'candidates':
                file_name = ('SEDfix_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_'
                             'inclusiveGCclass_Jun21_phangshst_candidates_bcw_v1p2_IR4.fits' %
                             helper_func.FileTools.target_names_no_zeros(target))
            else:
                file_name = (('SEDfix_PHANGS_IR4_%s_NewModelsNBHaUnionHaFLAG91pc_inclusiveGCcc_'
                             'inclusiveGCclass_Jun21_phangs_hst_v1p2_%s_%s.fits') %
                             (helper_func.FileTools.target_names_no_zeros(target), classify, cl_class))
        elif table_number == 2:
            if cl_class == 'candidates':
                file_name = ('SEDfix_%s_NewModelsHSTHaUnionHaFLAG11pc_inclusiveGCcc_'
                             'inclusiveGCclass_Jun21_phangshst_candidates_bcw_v1p2_IR4.fits' %
                             helper_func.FileTools.target_names_no_zeros(target))
            else:
                file_name = (('SEDfix_PHANGS_IR4_%s_NewModelsHSTHaUnionHaFLAG11pc_'
                              'inclusiveGCcc_inclusiveGCclass_Jun21_phangs_hst_v1p2_%s_%s.fits') %
                             (helper_func.FileTools.target_names_no_zeros(target), classify, cl_class))
        elif table_number == 3:
            file_name = '%s_HST_Ha_all_clusters_results.fits' % helper_func.FileTools.target_names_no_zeros(target)
        elif table_number == 4:
            file_name = ('%s_HST_Ha_nircam_all_clusters_results.fits' %
                         helper_func.FileTools.target_names_no_zeros(target))
        else:
            raise KeyError('Table number must be in 1,2,3,4')

        file_path_ir = path2table / file_name

        # check if file exists
        if not os.path.isfile(file_path_ir):
            print(file_path_ir, ' not found ! ')
            raise FileNotFoundError('there is no HST cluster catalog for the target ', target,
                                    ' make sure that the file ', file_path_ir, ' exists.')
        # open table and get column names
        table_hdu = fits.open(file_path_ir)
        table_ir = table_hdu[1].data
        table_hdu.close()

        return table_ir

    def get_artifact_cat(self, target):
        """
        Function to get artifact catalog name
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified.

        Returns
        -------
        table_artifact : ``astropy.io.fits.fitsrec.FITS_rec`` or None
            artifact table
        """
        # make sure the zero is removed after NGC for some objects
        if (target[:3] == 'ngc') & (target[3] == '0'):
            target_str = target[:3] + target[4:]
        else:
            target_str = target

        # get artifact table file path
        file_path_artifact = self.path2artifact / Path('%s_artifacts.fits' % target_str)
        # check if file exists
        if (not os.path.isfile(file_path_artifact)) & self.artifact_rais_file_not_found_flag:
            print(file_path_artifact, ' not found ! ')
            raise FileNotFoundError('there is no artigfact table for the target ', target,
                                    ' make sure that the file ', file_path_artifact, ' exists.')
        elif not os.path.isfile(file_path_artifact):
            print('No artifact table found for ', target)
            return None
        else:

            # open table and get column names
            table_hdu = fits.open(file_path_artifact)
            table_artifact = table_hdu[1].data
            table_hdu.close()

            return table_artifact

    def get_data_release_table_name(self, target, classify, cl_class, table_type='obs', sed_type=None):
        """
        Function to create final catalog name
        Parameters
        ----------
        target : str
            name of PHANGS-HST target. Must be in self.phangs_hst_target_list
        classify : str
            classification `human` or `ml`
        cl_class : str
            class group specification either `class12` or `class3`
        table_type : str
            to distinguish between obs and sed table
        sed_type : str
        Returns
        -------
        table_name : str
            table name
        """

        cat_str = 'hlsp_phangs-cat_hst_'

        # get instruments

        band_list = helper_func.BandTools.get_hst_obs_band_list(target=target, phangs_hst_version='phangs_hst_1')
        instrument_list = []
        for band in band_list:
            instrument_list.append(helper_func.BandTools.get_hst_instrument(target=target, band=band))

        instruments = ''
        if 'acs' in instrument_list:
            instruments += 'acs'
            if 'uvis' in instrument_list:
                instruments += '-uvis'
        else:
            instruments += 'uvis'

        cat_str += instruments + '_'
        cat_str += target.lower() + '_'
        cat_str += 'multi' + '_'
        cat_str += self.cat_ver + '_'
        if table_type == 'obs':
            cat_str += 'obs' + '-'
        if table_type == 'sed':
            cat_str += 'sed' + '-'
        if table_type == 'cand':
            cat_str += 'obs-sed' + '-'
        if classify == 'human':
            cat_str += 'human' + '-'
        if classify == 'ml':
            cat_str += 'machine' + '-'
        if cl_class == 'class12':
            cat_str += 'cluster-class12'
        if cl_class == 'class3':
            cat_str += 'compact-association-class3'
        if cl_class is None:
            cat_str += 'candidates'
        # add catalog_description
        if sed_type is not None:
            cat_str += '_' + sed_type
        cat_str += '.fits'

        return cat_str

    def create_obs_table(self, target, table, cand_table, classify):
        """
        Function to convert an internal data release table into a final data release table
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table
        cand_table : type ``astropy.io.fits.fitsrec.FITS_rec``
            candiate table for index crossmatch fits table
        classify : str
            classification human or ml
        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """
        # get column names for the catalog
        column_name_list = self.get_obs_table_column_list(target)
        # create table
        obs_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            # flux correction
            # if target in ['ngc1512', 'ngc1510']:
            #     if (self.cat_info[col_name]['col_name'][-3:] == 'mJy') | (self.cat_info[col_name]['col_name'][-7:] == 'mJy_ERR'):
            #         # now check which entrances have a detection
            #         mask_content = column_content != -9999.0
            #         column_content[mask_content] *= self.ngc1512_app_corr_offset_flux
            #     if self.cat_info[col_name]['col_name'][-4:] == 'VEGA':
            #         mask_content = column_content != -9999.0
            #         column_content[mask_content] += self.ngc1512_app_corr_offset_mag

            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment'])

            if obs_table is None:
                obs_table = column
            else:
                obs_table = hstack([obs_table, column])
        # now add cluster and all source id
        data_cluster_id = np.ones(len(obs_table), dtype=int) * -999
        data_allsource_id = np.ones(len(obs_table), dtype=int) * -999
        for obj_index in range(len(obs_table)):
            obj_mask = ((cand_table['PHANGS_X'] == obs_table['PHANGS_X'][obj_index]) &
                        (cand_table['PHANGS_Y'] == obs_table['PHANGS_Y'][obj_index]))
            data_cluster_id[obj_index] = cand_table['ID_PHANGS_CLUSTER'][obj_mask]
            data_allsource_id[obj_index] = cand_table['ID_PHANGS_ALLSOURCES'][obj_mask]

        cluster_id_column = Column(data=data_cluster_id,
                                   name=self.cat_info['id_phangs_cluster']['col_name'],
                                   dtype=data_cluster_id.dtype,
                                   description=self.cat_info['id_phangs_cluster']['doc_comment'])
        allsource_id_column = Column(data=data_allsource_id,
                                     name=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['col_name'],
                                     dtype=data_allsource_id.dtype,
                                     description=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['doc_comment'])

        obs_table.add_column(cluster_id_column, index=1)
        obs_table.add_column(allsource_id_column, index=3)

        # add region column:
        young_mask = table['SEDfix_age'] < 10
        vi_color = table['PHANGS_F555W_VEGA_TOT'] - table['PHANGS_F814W_VEGA_TOT']
        color_u = table['PHANGS_F336W_VEGA_TOT']
        if 'F438W' in self.phangs_hst_obs_band_dict[target]['wfc3_uvis_observed_bands']:
            color_b = table['PHANGS_F438W_VEGA_TOT']
        else:
            color_b = table['PHANGS_F435W_VEGA_TOT']
        ub_color = color_u - color_b

        vi_hull_ycl, ub_hull_ycl = self.load_hull(region_type='ycl', classify=classify, y_color='ub')
        vi_hull_map, ub_hull_map = self.load_hull(region_type='map', classify=classify, y_color='ub')
        vi_hull_ogcc, ub_hull_ogcc = self.load_hull(region_type='ogcc', classify=classify, y_color='ub')

        hull_ycl = ConvexHull(np.array([vi_hull_ycl, ub_hull_ycl]).T)
        hull_map = ConvexHull(np.array([vi_hull_map, ub_hull_map]).T)
        hull_ogcc = ConvexHull(np.array([vi_hull_ogcc, ub_hull_ogcc]).T)

        in_hull_ycl = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ycl)
        in_hull_map = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_map)
        in_hull_ogcc = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ogcc)

        mask_ycl = in_hull_ycl * np.invert(in_hull_map) + young_mask * (in_hull_map * in_hull_ogcc)
        mask_map = in_hull_map * np.invert(young_mask)
        mask_ogcc = in_hull_ogcc * np.invert(young_mask)

        column_data_region = np.array(['outside'] * len(table['SEDfix_age']), dtype=str)

        column_data_region[mask_ycl] = 'YCL'
        column_data_region[mask_map] = 'MAP'
        column_data_region[mask_ogcc] = 'OGCC'

        region_column = Column(data=column_data_region,
                               name=self.cat_info['cc_class']['col_name'],
                               dtype=str,
                               description=self.cat_info['cc_class']['doc_comment'])
        obs_table.add_column(region_column, index=-1)

        return obs_table

    def create_sed_table(self, target, table, cand_table):
        """
        Function to convert an internal data release table into a final data release table
        ----------
        target : str
            target for which the IR catalog should be identified. Must be in self.phangs_hst_target_list
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table
        cand_table : type ``astropy.io.fits.fitsrec.FITS_rec``
            candiate table for index crossmatch fits table
        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """
        # get column names for the catalog
        column_name_list = self.tab2_columns
        # create table
        sed_table = None
        for col_name in column_name_list:
            column_content = table[col_name]
            # mass correction
            # if target in ['ngc1512', 'ngc1510']:
            #     if (self.cat_info[col_name]['col_name'] in ['PHANGS_MASS_MINCHISQ', 'PHANGS_MASS_MINCHISQ_ERR']) | ('mass' in self.cat_info[col_name]['col_name']):
            #         # now check which entrances have a detection
            #         print(self.cat_info[col_name]['col_name'])
            #         mask_content = (column_content != -9999.0) & (column_content != -999.0)
            #         column_content[mask_content] *= self.ngc1512_app_corr_offset_flux

            if self.cat_info[col_name]['unit'] is not None:
                column_content *= self.cat_info[col_name]['unit']
            column = Column(data=column_content,
                            name=self.cat_info[col_name]['col_name'],
                            dtype=column_content.dtype,
                            description=self.cat_info[col_name]['doc_comment']
                            )

            if sed_table is None:
                sed_table = column
            else:
                sed_table = hstack([sed_table, column])

        # now add cluster and all source id
        data_cluster_id = np.ones(len(sed_table), dtype=int) * -999
        data_allsource_id = np.ones(len(sed_table), dtype=int) * -999
        for obj_index in range(len(sed_table)):
            obj_mask = ((cand_table['PHANGS_X'] == sed_table['PHANGS_X'][obj_index]) &
                        (cand_table['PHANGS_Y'] == sed_table['PHANGS_Y'][obj_index]))
            data_cluster_id[obj_index] = cand_table['ID_PHANGS_CLUSTER'][obj_mask]
            data_allsource_id[obj_index] = cand_table['ID_PHANGS_ALLSOURCES'][obj_mask]

        cluster_id_column = Column(data=data_cluster_id,
                                   name=self.cat_info['id_phangs_cluster']['col_name'],
                                   dtype=data_cluster_id.dtype,
                                   description=self.cat_info['id_phangs_cluster']['doc_comment'])
        allsource_id_column = Column(data=data_allsource_id,
                                     name=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['col_name'],
                                     dtype=data_allsource_id.dtype,
                                     description=self.cat_info['ID_PHANGS_ALLSOURCES_v1p2']['doc_comment'])

        sed_table.add_column(cluster_id_column, index=1)
        sed_table.add_column(allsource_id_column, index=3)

        return sed_table

    def create_cand_table(self, target, table):
        """
        Function to convert an internal data release table into a final data release table
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified.
        table : type ``astropy.io.fits.fitsrec.FITS_rec``
            input fits table

        Returns
        -------
        table_ir : ``astropy.table.Table``
            Final data release table for one object
        """

        # get column names which will enter into the catalog
        col_table = ascii.read(self.path2col_name_tab_1, format='csv')
        # use only those columns which will be in
        col_table = col_table[np.array(col_table['cand_tab'], dtype=bool)]

        # get the band name of the u-band
        band_list = helper_func.BandTools.get_hst_obs_band_list(target=target, phangs_hst_version='phangs_hst_1')
        if 'F435W' in band_list:
            not_u_band = 'F438W'
        elif 'F438W' in band_list:
            not_u_band = 'F435W'
        else:
            raise RuntimeError('Something went wrong. For the galaxy ', target, ' There is no U-band available! It must be F435W or F438W ')

        cand_table = None
        for raw_col_name, final_col_name, dtype, astropy_unit_str, doc_description in zip(
                col_table['Provided_col_name'], col_table['final_col_name'], col_table['dtype'],
                col_table['astropy_unit_str'], col_table['doc_description']):
            # kick out the wrong U-band name
            if not_u_band in raw_col_name:
                continue

            column_content = table[raw_col_name]
            # add astropy unit if possible
            if isinstance(astropy_unit_str, str):
                column_content *= u.Unit(astropy_unit_str)

            column = Column(data=column_content, name=final_col_name, dtype=dtype.astype('str'),
                            description=str(doc_description))
            if cand_table is None:
                cand_table = column
            else:
                cand_table = hstack([cand_table, column])

        return cand_table

    def load_second_classify_table(self):
        """
        Function to load the secondary classification table into the attributes
        Parameters
        ----------
        ...
        Returns
        -------
        None
        """
        self.questionable_artefact_table = ascii.read(self.path2questionable_artifacts)

    def get_second_classify_cat(self, target):
        """
        Function to get artifact catalog that was reclassified because the first artefact classification was
        questionable.
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified.

        Returns
        -------
        table_artifact : ``astropy.io.fits.fitsrec.FITS_rec`` or None
            artifact table
        """
        # make sure the zero is removed after NGC for some objects

        if self.questionable_artefact_table is None:
            self.load_second_classify_table()

        if (target[:3] == 'ngc') & (target[3] == '0'):
            target_str = target[:3] + target[4:]
        else:
            target_str = target

        mask_target = self.questionable_artefact_table['GALAXY'] == target_str

        return self.questionable_artefact_table[mask_target]

    def load_diffraction_spike_masks(self, target, target_str):
        """
        Function to get masks for diffraction spikes, produced by Kirsten
        Parameters
        ----------
        target : str
            target for which the IR catalog should be identified.

        Returns
        -------
        diffraction_spike_mask : ndarray
        """

        if (target_str[:3] == 'ngc') & (target_str[3] == '0'):
            target_str = target[:3] + target[4:]

        diffraction_mask = None
        diffraction_wcs = None
        # we just sum up all maps (0 nothing and 1 indicates diffraction spike)
        for band in helper_func.BandTools.get_hst_obs_band_list(target=target, phangs_hst_version='phangs_hst_1'):
            instrument = helper_func.BandTools.get_hst_instrument(target=target, band=band)

            diff_spike_mask_file_name = (self.path2diffraction_spike_masks + '/' '%s_%s_%s_mask.fits' %
                                         (target_str, band.lower(), instrument))
            if os.path.isfile(diff_spike_mask_file_name):
                hdu = fits.open(diff_spike_mask_file_name)
                if diffraction_mask is None:
                    diffraction_mask = hdu[0].data
                    diffraction_wcs = WCS(hdu[0].header)
                else:
                    diffraction_mask += hdu[0].data

        return diffraction_mask, diffraction_wcs

    def load_hull(self, region_type, classify, y_color='ub'):
        if region_type == 'ycl':
            region_str = 'young'
            class_number = 3
        elif region_type == 'map':
            region_str = 'mid'
            class_number = 1
        elif region_type == 'ogcc':
            region_str = 'ogc'
            class_number = 1
        else:
            raise KeyError('region string not understood')

        if classify == 'human':
            classify_str = 'hum'
        elif classify == 'ml':
            classify_str = 'ml'
        else:
            raise KeyError('classify string not understood')

        x_hull = np.load('ccd_region_hull/vi_hull_%s_%svi_%s_%i.npy' %
                         (region_str, y_color, classify_str, class_number))
        y_hull = np.load('ccd_region_hull/%s_hull_%s_%svi_%s_%i.npy' %
                         (y_color, region_str, y_color, classify_str, class_number))

        return x_hull, y_hull

    def add_ccd_region_col(self, target, candidate_table):
        """
        Function to add CCD region classifications
        """
        # add column for CCD regions (See Maschmann+24 for more details)
        # add region column:
        young_mask = candidate_table['PHANGS_SED_AGE'] < 10
        vi_color = candidate_table['PHANGS_F555W_VEGA'] - candidate_table['PHANGS_F814W_VEGA']
        color_u = candidate_table['PHANGS_F336W_VEGA']

        hst_band_list = helper_func.BandTools.get_hst_obs_band_list(target=target, phangs_hst_version='phangs_hst_1')
        if 'F438W' in hst_band_list:
            color_b = candidate_table['PHANGS_F438W_VEGA']
        else:
            color_b = candidate_table['PHANGS_F435W_VEGA']
        ub_color = color_u - color_b

        # load hulls
        vi_hull_ycl_hum, ub_hull_ycl_hum = self.load_hull(region_type='ycl', classify='human', y_color='ub')
        vi_hull_map_hum, ub_hull_map_hum = self.load_hull(region_type='map', classify='human', y_color='ub')
        vi_hull_ogcc_hum, ub_hull_ogcc_hum = self.load_hull(region_type='ogcc', classify='human', y_color='ub')
        vi_hull_ycl_ml, ub_hull_ycl_ml = self.load_hull(region_type='ycl', classify='ml', y_color='ub')
        vi_hull_map_ml, ub_hull_map_ml = self.load_hull(region_type='map', classify='ml', y_color='ub')
        vi_hull_ogcc_ml, ub_hull_ogcc_ml = self.load_hull(region_type='ogcc', classify='ml', y_color='ub')

        hull_ycl_hum = ConvexHull(np.array([vi_hull_ycl_hum, ub_hull_ycl_hum]).T)
        hull_map_hum = ConvexHull(np.array([vi_hull_map_hum, ub_hull_map_hum]).T)
        hull_ogcc_hum = ConvexHull(np.array([vi_hull_ogcc_hum, ub_hull_ogcc_hum]).T)
        hull_ycl_ml = ConvexHull(np.array([vi_hull_ycl_ml, ub_hull_ycl_ml]).T)
        hull_map_ml = ConvexHull(np.array([vi_hull_map_ml, ub_hull_map_ml]).T)
        hull_ogcc_ml = ConvexHull(np.array([vi_hull_ogcc_ml, ub_hull_ogcc_ml]).T)

        in_hull_ycl_hum = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ycl_hum)
        in_hull_map_hum = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_map_hum)
        in_hull_ogcc_hum = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ogcc_hum)
        in_hull_ycl_ml = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ycl_ml)
        in_hull_map_ml = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_map_ml)
        in_hull_ogcc_ml = helper_func.points_in_hull(np.array([vi_color, ub_color]).T, hull_ogcc_ml)

        mask_ycl_hum = in_hull_ycl_hum * np.invert(in_hull_map_hum) + young_mask * (in_hull_map_hum * in_hull_ogcc_hum)
        mask_map_hum = in_hull_map_hum * np.invert(young_mask)
        mask_ogcc_hum = in_hull_ogcc_hum * np.invert(young_mask)
        mask_ycl_ml = in_hull_ycl_ml * np.invert(in_hull_map_ml) + young_mask * (in_hull_map_ml * in_hull_ogcc_ml)
        mask_map_ml = in_hull_map_ml * np.invert(young_mask)
        mask_ogcc_ml = in_hull_ogcc_ml * np.invert(young_mask)

        column_data_region_hum = np.array(['outside'] * len(candidate_table['PHANGS_SED_AGE']), dtype=str)
        column_data_region_ml = np.array(['outside'] * len(candidate_table['PHANGS_SED_AGE']), dtype=str)

        column_data_region_hum[mask_ycl_hum] = 'YCL'
        column_data_region_hum[mask_map_hum] = 'MAP'
        column_data_region_hum[mask_ogcc_hum] = 'OGCC'
        column_data_region_ml[mask_ycl_ml] = 'YCL'
        column_data_region_ml[mask_map_ml] = 'MAP'
        column_data_region_ml[mask_ogcc_ml] = 'OGCC'

        region_column_hum = Column(data=column_data_region_hum, name='CC_CLASS_HUMAN', dtype=str,
                                   description='Flag to identify in which region on the color-color diagram the object '
                                               'was associated with. Values are \"ycl\" (young cluster locus), \"map\" '
                                               '(middle aged plume) \"ogcc\" (old globular cluster clump) or '
                                               '\"outside\" (outside the main regions and therefore not classified). '
                                               'A detailed description is found in Section 4.4 of Maschmann+24')
        region_column_ml = Column(data=column_data_region_ml, name='CC_CLASS_ML', dtype=str,
                                  description='Flag to identify in which region on the color-color diagram the object '
                                              'was associated with. Values are \"ycl\" (young cluster locus), \"map\" '
                                              '(middle aged plume) \"ogcc\" (old globular cluster clump) or '
                                              '\"outside\" (outside the main regions and therefore not classified). '
                                              'A detailed description is found in Section 4.4 of Maschmann+24')
        index_last_classification_hum = np.where(np.array(candidate_table.colnames) == 'PHANGS_CI')
        candidate_table.add_column(region_column_hum, index=index_last_classification_hum[0][0]+1)

        index_last_classification_ml = np.where(np.array(candidate_table.colnames) == 'CC_CLASS_HUMAN')
        candidate_table.add_column(region_column_ml, index=index_last_classification_ml[0][0]+1)