"""
A module to download the truth parameters of astrophysical objects on fatboy and
build light curves from them.
"""

from __future__ import absolute_import, division, print_function
import numpy as np
import pandas as pd
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.supernovae import SNObject
import pymssql
import os
import json
from collections import OrderedDict
from lsst.utils import getPackageDir
from lsst.daf.persistence import DbAuth
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import Variability as var
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator as obsGen


__all__ = ['RefLightCurves_SL']


class RefLightCurves_SL(object):
    """
    Class to connect to the database tables on fatboy and obtain truth
    light curves

    Parameters
    ----------
    tableName : string, mandatory
        case insensitive string name of table on database to connect to
        for model parameters of astrophysical objects
    idCol : string, optional, defaults to 'snid'
        column name of Index on the table
    observations : `pd.DataFrame`, optional, defaults to None
        if None, some information may need to be supplied in using certain
        methods like lightCurve. The dataframe must have the following columns
        'index' : obsHistID, 'expMJD': MJD of the time of exposure,
        'fiveSigmaDepth': magnitude of a point source for which the signal to
        noise ratio in that exposure would be expected to be 5.0. Such
        information would be accessible from a OpSim output.
    bandPassDict : instance of `lsst.sims.bandPassDict`, optional,
        defaults to None
    columns : tuple of strings, optional, defaults to values for SN
        tuple of strings that completely specify the truth values for
        the astrophysical object
    dbConnection : `pymssql.connection` instance, mandatory
        connection to the database where the relevant tables of catsim
        objects are stored
    dbCursor : `pymssql.connection.cursor` instnce, optional, defaults to None
        cursor to the DataBase Connection. If None, a new cursor is obtained
        from self.dbConnection
    idSequence : sequence of one dimension, optional, defaults to None
        sequence of unique ids in the catsim universe indexing the
        astrophysical objects in the database.
    dbHostName : string, optional, defaults to None
        force the class to use this hostname. If not provided, the class will
        set this to localhost, which is the desired hostname when using an ssh
        tunnel. This parameter is useful when working from whitelisted
        computers

    Attributes
    ----------
    columns :  The columns that will be obtained from the table on the catsim
        database
    idCol : The column name on the table of the astrophysical object which
        indexes the objects
    tableName : Name of the table on the catsim database containing the
        parameters of the astrophysical object
    objectTypeID : a unique integer associated with each class of astrophysical
        objects in catsim.
    dbHostName : a name for hostname of the database
    bandPassDict : bandpasses for the observations
    observations : a set of observations incorporating information about time
        of observaations and the fiveSigma Depths of the observations.

    Methods
    -------

    Examples
    --------
    >>> reflc = RefLightCurves(idSequence=(6001163623700, 6000324908000),
                               tableName='TwinkSN',
                               dbConnection=DBConnection,
                               dbCursor=db) # doctest: +SKIP
    """
    def __init__(self,
                 tableName,
                 objectTypeID=28,
                 idCol='galtileid',
                 columns=('galtileid', 'redshift', 'raJ2000', 'decJ2000',
                          'varParamStr'),
                 opsim_database=None,
                 observations=None,
                 bandPassDict=None,
                 dbConnection=None,
                 dbCursor=None,
                 dbHostName=None,
                 idSequence=None,
                 t0_mjd=None):

        self.columns = columns
        self._dbConnection = dbConnection
        self._dbCursor = dbCursor
        self._idSequence = idSequence
        self.columns = columns
        self.opsim_database = opsim_database
        self.idCol = idCol
        self.tableName = tableName
        self._idvals = None
        self.objectID = objectTypeID
        self.dbHostName = dbHostName
        self.bandPassDict = bandPassDict
        self.observations = observations
        self.t0_mjd = t0_mjd

    @property
    def dbConnection(self):
        """
        The pymssql connection to the catsim database used to query refrence
        objects
        """
        if self._dbConnection is None:
            config = bcm.BaseCatalogConfig()
            config.load(os.path.join(getPackageDir("sims_catUtils"), "config",
                                     "db.py"))

            username = DbAuth.username(config.host, config.port)
            password = DbAuth.password(config.host, config.port)
            hostname = config.host
            if self.dbHostName is not None:
                hostname = self.dbHostName
            DBConnection = pymssql.connect(user=username,
                                           password=password,
                                           host=hostname,
                                           database=config.database,
                                           port=config.port)
            return DBConnection
        else:
            return self._dbConnection

    @property
    def dbCursor(self):
        """
        Cursor to the catsim database connection. This is not reset if one
        exists.
        """
        if self._dbCursor is None:
            self._dbCursor = self.dbConnection.cursor()
        return self._dbCursor

    def _generateObsMetaData(self, ra, dec):
        obs_gen = obsGen(self.opsim_database)
        obs_metadata = obs_gen.getObservationMetaData(fieldRA=(ra-0.1, ra+0.1),
                                                      fieldDec=(dec-0.1, dec+0.1),
                                                      limit=1)
        return obs_metadata

    def _query_columns(self, obs_metadata):
        gto = bcm.GalaxyAgnObj()
        results = gto.query_columns(colnames=list(self.columns),
                                    obs_metadata=obs_metadata,
                                    chunk_size=None)
        ret_results = [tuple(xx) for xx in results][0]

        return ret_results

    def get_sl_params(self, ra, dec):
        obs_metadata = self._generateObsMetaData(ra, dec)
        results = self._query_columns(obs_metadata[0])

        results_dict = OrderedDict()
        for col_num, col_name in enumerate(self.columns):
            col_list = []
            if col_name == 'varParamStr':
                var_dict = OrderedDict()
                var_dict['t0_mjd'] = []
                var_dict['seed'] = []
                var_dict['agn_tau'] = []
                var_dict['agn_sfu'] = []
                var_dict['agn_sfg'] = []
                var_dict['agn_sfr'] = []
                var_dict['agn_sfi'] = []
                var_dict['agn_sfz'] = []
                var_dict['agn_sfy'] = []
            for res_line in results:
                if col_name == 'varParamStr':
                    varInfo = json.loads(res_line[col_num])
                    varParams = varInfo['pars']
                    for k, v in varParams.items():
                        var_dict[k].append(v)
                    if self.t0_mjd is not None:
                        for idx in xrange(len(var_dict['t0_mjd'])):
                            var_dict['t0_mjd'][idx] = self.t0_mjd
                else:
                    if (col_name == ('raJ2000') or col_name == ('decJ2000')):
                        col_list.append(np.degrees(res_line[col_num]))
                    else:
                        col_list.append(res_line[col_num])
            if col_name != 'varParamStr':
                results_dict[col_name] = col_list

        df_output = pd.DataFrame.from_dict(results_dict)
        df_var = pd.DataFrame.from_dict(var_dict)

        return pd.concat([df_output, df_var], axis=1)

    @staticmethod
    def uniqueIDtoTableId(uniqueID, objTypeID, nshift=10):
        """
        Given a sequence of catsim uniqueIDs, convert it to a numpy
        array of IDs in the table of the object (called refIDCol) using
        objTypeID.

        Parameters
        ----------
        uniqueID : 1D sequence of unique IDs as found in catsim/phosim Instance
            catalogs.
        objTypeID : A unique ID assigned to each class of object in the catsim
            database.
        nshift : integer, optional, defaults to 10
            Number of bit shifts, exactly the same as in catsim.

        Returns
        -------
        `numpy.ndarray` of IDs indexing the table of the particular object.

        .. note:: This is an inverse of the catsim function
            `lsst.sims.catalogs_measures.Instance.get_uniqueId`. Later on I
            hope this code will be moved to a similar location.
        """
        id = np.asarray(uniqueID) - objTypeID
        return np.right_shift(id, nshift)

    def agnLightCurve(self, df_row, mjdList):

        agn_lc = var()
        lightCurveDicts = []
        for mjd in mjdList:
            dmag = agn_lc.applyAgn(df_row, mjd)
            dmag['mjd'] = mjd
            lightCurveDicts.append(dmag)
        agn_lc_df = pd.DataFrame(lightCurveDicts)

        return agn_lc_df


    def lightCurve(self, idValue,
                   bandName=None,
                   observations=None,
                   bandPassDict=None,
                   format='dataframe',
                   keys=['expMJD', 'filter', 'fiveSigmaDepth'],
                   photParams=None):
        """
        return the light curve of the object for observations as a
        `pandas.Dataframe`.
        Parameters
        ----------
        idValue : integer, mandatory
            index of astrophysical object
        bandName : string, optional, defaults to None
            key for bandpassDict, so that bandpassDict[key] is an instance
            of `lsst.sims.photUtils.Bandpass`. If provided, only the light
            curve in band band is returned
        observations : optional, allowed formats decided by format variable,
            defaults to None
            maximal set observations corresponding to which the light curve is
            obtained. If None, observations default to self.observaions. If
            both observations and self.observations are None, an exception
            will be raised.
        bandPassDict : instance of `lsst.sims.photUtils.BandPassDict`, optional
            defaults to None
            dictionary of total (system + atmospheric) bandpasses as a
            dictionary of bands. If None, defaults to self.bandPassDict. If
            self.bandPassDict is None, an exception is raised
        format: string, optional, defaults to dataframe
            format in which observations are available
        keys: aliases for columns time (in MJD), bandpass Name, fivesigmaDepth
            (mags)
        photParams: instance of `sims.photUtils.PhotParams`, optional,
            defaults to None
            Describes the observing conditions and telescope for the
            observation. The default value of None instantiates this for the
            LSST site and telescope.

        Returns
        -------
        `pd.DataFrame` containing the following columns at the minimum
        ['time', 'flux', 'fluxerr', 'band']. time is in modified Julian Days,
        flux and fluxerrr are in `maggies`. If a set of indexes are provided
        for the m5Values (eg. as indexes in a dataFrame)
        """
        format = format.lower()

        if format not in ['dataframe']:
            raise NotImplementedError(
                "Unavailable input format {0}".format(format))

        sn = self.astro_object(idValue)
        timeMin = sn.mintime()
        timeMax = sn.maxtime()

        if bandPassDict is None:
            bandPassDict = self.bandPassDict
        if bandPassDict is None:
            raise ValueError('The method parameter bandPassDict, and '
                             'the attribute bandPassDict cannot simultaneously'
                             ' be None\n')
        if observations is None:
            observations = self.observations
        if observations is None:
            raise ValueError('The method parameter observations, and '
                             'the attribute observations cannot simultaneously'
                             ' be None\n')
        df = observations.query('expMJD < @timeMax and '
                                'expMJD > @timeMin').copy()
        times = []
        bands = []
        fluxs = []
        fluxerrs = []
        m5vals = []
        for ind in df.index.values:
            time = df.ix[ind, 'expMJD']
            band = df.ix[ind, 'filter']
            flux = sn.catsimBandFlux(time=time,
                                     bandpassobject=bandPassDict[band])
            m5val = df.ix[ind, 'fiveSigmaDepth']
            fluxerr = sn.catsimBandFluxError(time=time,
                                             bandpassobject=bandPassDict[band],
                                             m5=m5val, fluxinMaggies=flux)
            times.append(time)
            bands.append(band)
            m5vals.append(m5val)
            fluxs.append(flux)
            fluxerrs.append(fluxerr)

        mydict = dict()
        mydict['time'] = times
        mydict['flux'] = fluxs
        mydict['fluxerr'] = fluxerrs
        mydict['band'] = bands
        mydict['m5'] = m5vals
        output = pd.DataFrame(mydict, index=df.index)
        if bandName is not None:
            output = output.query('band == @bandName')
        return output
