from __future__ import absolute_import, division, print_function
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import (DefaultPhoSimHeaderMap,
                                                         PhoSimCatalogPoint)
from lsst.sims.photUtils import BandpassDict, SedList

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
               ('sedFilename', 'sedfilename', unicode, 40)]

class TrueStars(object):

    def __init__(self, dbConn, opsimDB_filename):

        self.dbConn = dbConn
        self.opsimDB = opsimDB_filename
        # Set up OpSim database (from Twinkles/bin/generatePhosimInput.py)
        engine = create_engine('sqlite:///' + self.opsimDB)
        self.obs_gen = ObservationMetaDataGenerator(database=self.opsimDB,
                                                    driver='sqlite')

    def get_true_stars(self, for_obsHistIds=None):

        if for_obsHistIds is None:
            raise TypeError("Please specify a visit using 'for_visit='.")

        obs_metadata_list = []
        for obsHistID in for_obsHistIds:
            obs_metadata_list.append(self.obs_gen.getObservationMetaData(
                                                    obsHistID=obsHistID,
                                                    fieldRA=(53, 54),
                                                    fieldDec=(-29, -27),
                                                    boundLength=0.3)[0])

        star_df = pd.DataFrame(columns = ['uniqueId', 'filter',
                                          'true_flux', 'obsHistId'])
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        bp_indices = {}
        for bp in list(enumerate(bp_dict.keys())):
            bp_indices[bp[1]] = bp[0]

        column_names = None
        self.obs_md = obs_metadata_list[0]

        for obs_metadata in obs_metadata_list:
            star_cat = PhoSimCatalogPoint(self.dbConn, obs_metadata=obs_metadata)
            if column_names is None:
                column_names = [x for x in star_cat.iter_column_names()]
            star_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
            self.star_cat = star_cat
            chunk_data = []
            i = 0
            for line in self.star_cat.iter_catalog():
                chunk_data.append(line)
            chunk_data = pd.DataFrame(chunk_data, columns=column_names)

            sed_list = SedList(chunk_data['sedFilepath'],
                               chunk_data['phoSimMagNorm'],
                               specMap = None,
                               galacticAvList=chunk_data['galacticAv'])
            mag_list = bp_dict.fluxListForSedList(sed_list,
                          indices=[bp_indices[obs_metadata.OpsimMetaData['filter']]])
            mag_list = np.array(mag_list)[:,bp_indices[obs_metadata.OpsimMetaData['filter']]]
            visit_df = pd.DataFrame(np.array([chunk_data['uniqueId'],
                                    [obs_metadata.OpsimMetaData['filter']]*len(chunk_data),
                                    mag_list,
                                    [obs_metadata.OpsimMetaData['obsHistID']]*len(chunk_data)]).T,
                                    columns = ['uniqueId', 'filter',
                                               'true_flux', 'obsHistId'])
            star_df = star_df.append(visit_df, ignore_index=True)

        self.star_df = star_df
