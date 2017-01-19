from __future__ import absolute_import, division, print_function
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.photUtils import BandpassDict, SedList
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from desc.monitor.TruthCatalogDefs import TruthCatalogPoint

__all__ = ["StarCacheDBObj", "TrueStars"]

class StarCacheDBObj(CatalogDBObject):
    tableid = 'star_cache_table'
    host = None
    port = None
    driver = 'sqlite'
    objectTypeId = 4
    idColKey = 'simobjid'
    raColName = 'ra'
    decColName = 'decl'

    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 256)]

class TrueStars(object):

    def __init__(self, dbConn, opsimDB_filename):

        self.dbConn = dbConn
        self.opsimDB = opsimDB_filename
        # Set up OpSim database (from Twinkles/bin/generatePhosimInput.py)
        engine = create_engine('sqlite:///' + self.opsimDB)
        self.obs_gen = ObservationMetaDataGenerator(database=self.opsimDB,
                                                    driver='sqlite')

    def get_true_stars(self, for_obsHistIds=None):

        """
        Get all the fluxes for stars in all visits in Twinkles.

        Can specify a subset of visits with for_obsHistIds.
        """

        if for_obsHistIds is None:
            survey_info = np.genfromtxt('../data/selectedVisits.csv',
                                    names=True, delimiter=',')
            for_obsHistIds = survey_info['obsHistID']

        obs_metadata_list = []
        visit_on = 0
        for obsHistID in for_obsHistIds:
            if visit_on % 100 == 0:
                print("Generated %i out of %i obs_metadata" %
                      (visit_on+1, len(for_obsHistIds)))
            visit_on += 1
            obs_metadata_list.append(self.obs_gen.getObservationMetaData(
                                                    obsHistID=obsHistID,
                                                    fieldRA=(53, 54),
                                                    fieldDec=(-29, -27),
                                                    boundLength=0.3)[0])

        star_df = pd.DataFrame(columns = ['uniqueId', 'ra', 'dec', 'filter',
                                          'true_flux', 'true_flux_error',
                                          'obsHistId'])
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        bp_indices = {}
        for bp in list(enumerate(bp_dict.keys())):
            bp_indices[bp[1]] = bp[0]

        column_names = None
        seds_loaded = False

        visit_on = 0
        for obs_metadata in obs_metadata_list:
            if visit_on % 100 == 0:
                print("Generated fluxes for %i out of %i visits" %
                      (visit_on+1, len(for_obsHistIds)))
            visit_on += 1
            star_cat = TruthCatalogPoint(self.dbConn, obs_metadata=obs_metadata,
                                         constraint='gmag > 11')

            if column_names is None:
                column_names = [x for x in star_cat.iter_column_names()]
            star_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
            chunk_data = []
            for line in star_cat.iter_catalog():
                chunk_data.append(line)
            chunk_data = pd.DataFrame(chunk_data, columns=column_names)

            #All SEDs will be the same since we are looking at the same point
            #in the sky and mag_norms will be the same for stars.
            if seds_loaded is False:
                sed_list = SedList(chunk_data['sedFilepath'],
                                   chunk_data['phoSimMagNorm'],
                                   specMap = None,
                                   galacticAvList=chunk_data['galacticAv'])
                seds_loaded = True

                mag_array = bp_dict.magArrayForSedList(sed_list)
                phot_params = PhotometricParameters()

            visit_filter = obs_metadata.OpsimMetaData['filter']
            flux_array = np.power(10,-0.4*(mag_array[visit_filter] - 22.5))
            snr, gamma = calcSNR_m5(mag_array[visit_filter],
                                    bp_dict[visit_filter],
                                    obs_metadata.OpsimMetaData['fiveSigmaDepth'],
                                    phot_params)
            flux_error = flux_array/snr

            visit_df = pd.DataFrame(np.array([chunk_data['uniqueId'],
                                    chunk_data['raJ2000'],
                                    chunk_data['decJ2000'],
                                    [visit_filter]*len(chunk_data),
                                    flux_array, flux_error,
                                    [obs_metadata.OpsimMetaData['obsHistID']]*len(chunk_data)]).T,
                                    columns = ['uniqueId', 'ra', 'dec', 'filter',
                                               'true_flux', 'true_flux_error',
                                               'obsHistId'])
            star_df = star_df.append(visit_df, ignore_index=True)

        self.star_df = star_df

    def write_to_db(self, filename, table_name='stars'):

        disk_engine = create_engine('sqlite:///%s' % filename)
        worker.star_df.to_sql('stars', disk_engine)
