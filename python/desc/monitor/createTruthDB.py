from __future__ import absolute_import, division, print_function
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.photUtils import BandpassDict, SedList
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from desc.monitor.truthCatalogDefs import TruthCatalogPoint

__all__ = ["StarCacheDBObj", "TrueStars"]


class StarCacheDBObj(CatalogDBObject):
    """
    CatalogDBObject for the stars in the simulation "truth" database.
    """
    tableid = 'star_cache_table'
    host = None
    port = None
    driver = 'sqlite'
    objectTypeId = 4
    idColKey = 'simobjid'
    raColName = 'ra'
    decColName = 'decl'

    columns = [('id', 'simobjid', int),
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
    """
    Gets the stars out of the simulation "truth" database for a specified set
    of visits.

    "True" in this case refers to the values that come from our LSST                                                                                            
    CATSIM database. This CATSIM database stores the LSST simulated universe                                                                                    
    model that we use to provide inputs to LSST simulations. It is important                                                                                    
    to note that this means that "true" does not refer to actual stars in                                                                                       
    sky, but to the known inputs to our simulations. More information                                                                                           
    on the LSST Simulations can be found here: bit.ly/lsst-sims-doc. 

    Note : RA, DEC values in "truth" catalogs are J2000 coordinates. Flux
    values in final output here are in nanomaggies.

    Parameters
    ----------
    dbConn : dbInterface instance
        This is a connection to a cached database of the simulation inputs.
    opsimDB_filename : str
        The location of the opsim database used with the simulation.
    """

    def __init__(self, dbConn, opsimDB_filename):

        self.dbConn = dbConn
        self.opsimDB = opsimDB_filename
        # Set up OpSim database (from Twinkles/bin/generatePhosimInput.py)
        self.obs_gen = ObservationMetaDataGenerator(database=self.opsimDB,
                                                    driver='sqlite')

    def get_true_stars(self, for_obsHistIds=None, catalog_constraints=None):
        """
        Get all the fluxes for stars in all visits specified.

        Parameters
        ----------
        for_obsHistIds : int or list of ints or None, default=None
            Can specify a subset of visits. If set to None will get it from the
            file in the repo located at monitor/data/selectedVisits.csv
        catalog_constraints : str or None, default=None
            Specify SQL constraints on the sims catalog used as
            the "truth" input.

        Returns
        ----------
        star_df : pandas dataframe
            Stores all the star information for the simulation inputs across
            the desired visits
        """

        if for_obsHistIds is None:
            survey_info = np.genfromtxt('../data/selectedVisits.csv',
                                        names=True, delimiter=',')
            for_obsHistIds = survey_info['obsHistID']
        else:
            for_obsHistIds = np.ravel(for_obsHistIds)

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

        star_df = pd.DataFrame(columns=['uniqueId', 'ra', 'dec', 'filter',
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
            star_cat = TruthCatalogPoint(self.dbConn,
                                         obs_metadata=obs_metadata,
                                         constraint=catalog_constraints)

            if column_names is None:
                column_names = [x for x in star_cat.iter_column_names()]
            star_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
            chunk_data = []
            for line in star_cat.iter_catalog():
                chunk_data.append(line)
            chunk_data = pd.DataFrame(chunk_data, columns=column_names)

            # All SEDs will be the same since we are looking at the same point
            # in the sky and mag_norms will be the same for stars.
            if seds_loaded is False:
                sed_list = SedList(chunk_data['sedFilepath'],
                                   chunk_data['phoSimMagNorm'],
                                   specMap=None,
                                   galacticAvList=chunk_data['galacticAv'])
                seds_loaded = True

                mag_array = bp_dict.magArrayForSedList(sed_list)
                flux_array = bp_dict.fluxArrayForSedList(sed_list)
                phot_params = PhotometricParameters()

            visit_filter = obs_metadata.OpsimMetaData['filter']
            # Get flux and convert to nanomaggies
            flux_array_visit = flux_array[visit_filter]/3.631e-06 
            five_sigma_depth = obs_metadata.OpsimMetaData['fiveSigmaDepth']
            snr, gamma = calcSNR_m5(mag_array[visit_filter],
                                    bp_dict[visit_filter],
                                    five_sigma_depth,
                                    phot_params)
            flux_error = flux_array_visit/snr

            obs_hist_id = obs_metadata.OpsimMetaData['obsHistID']
            visit_df = pd.DataFrame(np.array([chunk_data['uniqueId'],
                                    chunk_data['raJ2000'],
                                    chunk_data['decJ2000'],
                                    [visit_filter]*len(chunk_data),
                                    flux_array_visit, flux_error,
                                    [obs_hist_id]*len(chunk_data)]).T,
                                    columns=['uniqueId', 'ra', 'dec', 'filter',
                                             'true_flux', 'true_flux_error',
                                             'obsHistId'])
            star_df = star_df.append(visit_df, ignore_index=True)
            star_df['uniqueId'] = pd.to_numeric(star_df['uniqueId'])
            star_df['ra'] = pd.to_numeric(star_df['ra'])
            star_df['dec'] = pd.to_numeric(star_df['dec'])
            star_df['true_flux'] = pd.to_numeric(star_df['true_flux'])
            t_f_error = 'true_flux_error'
            star_df[t_f_error] = pd.to_numeric(star_df[t_f_error])
            star_df['obsHistId'] = pd.to_numeric(star_df['obsHistId'])

        return star_df

    def write_to_db(self, star_df, filename, table_name='stars', **kwargs):
        """
        Write self.star_df to a sqlite database.

        Parameters
        ----------
        star_df : pandas dataframe
            Stores all the star information for the simulation inputs across
            the desired visits        
        filename : str
            File name to use for the sqlite database.
        table_name : str, default='stars'
            Table name within the sqlite database for the star_df info.
        **kwargs
            Keyword arguments for the pandas `to_sql` function.
        """

        disk_engine = create_engine('sqlite:///%s' % filename)
        star_df.to_sql('stars', disk_engine, **kwargs)
