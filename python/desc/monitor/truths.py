"""
A class to download the truth parameters of objects and build light curves
from them.
"""

from __future__ import absolute_import, division, print_function

import os

import pandas as pd
import pymssql
from lsst.utils import getPackageDir
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.daf.persistence import DbAuth
from lsst.sims.catUtils.supernovae import SNObject


class RefLightCurves(object):
    """
    """
    def __init__(self,
                 idSequence,
                 dbCursor,
                 tableName,
                 dbConnection=None,
                 idCol='snid',
                 columns=('snid', 'redshift', 'snra', 'sndec', 't0', 'x0',
                          'x1', 'c')):
        """
        idSequence: sequence of integers, mandatory
            a sequence of SNIDs which index the SN on fatboy
        dbCursor:
            A curson to a DBConnection to a database containing the source SN
            properties.
        tableName:
        columns:
        """
        self.columns = columns
        self.dbConnection = dbConnection
        self.dbCursor = dbCursor
        self.idSequence = idSequence
        self.columns = columns
        self.idCol = idCol
        self.tableName = tableName

    def query(self, idValue=None):
        """
        Return the query statement to be used

        Parameters
        ----------
        idValue : integer, optional, defaults to None
            ID of the astrophysical object
        """
        query = """SELECT """
        query += ", ".join(xx.strip() for xx in self.columns)
        query += " FROM {} ".format(self.tableName)

        # if query is for a single idvalue, construct the query stmt
        if idValue is not None:
            query += "WHERE {0} = {1}".format(self.idCol, idValue)
        # if idValue is not supplied, but an idSequence is supplied
        elif self.idSequence is not None:
            query += "WHERE {0} in {1}".format(self.idCol, self.idSequence)
        # Else get the entire table, no WHERE clause
        else:
            pass
        return query

    def get_params(self, idValue=None):
        """
        return parameters of the objects with the desired columns as a
        `pandas.DataFrame`
        """
        df = pd.read_sql_query(self.query(idValue=idValue), self.dbConnection,
                               index_col=self.idCol, coerce_float=False)
        return df

    def astro_object(self, idValue, obsMetaDataList):
        """
        """
        df = self.get_params(idValue)
        sn = SNObject(ra=df.snra.values[0], dec=df.sndec.values[0])
        paramDict = dict()
        for param in ['t0', 'x0', 'x1', 'c']:
            paramDict[param] = df[param].values
        paramDict['z'] = df.redshift.values[0]
        sn.set(**paramDict)
        return sn

    def lightCurve(self, observations, format='dataframe',
                   keys=['expMJD', 'filter', 'fiveSigmaDepth'],
                   photParams=None):
        """
        return the light curve of the object for observations.

        Parameters
        ----------
        observations: mandatory, allowed formats decided by format variable
            observations corresponding to which the light curve is obtained
        format: string, optional, defaults to dataframe
        keys: aliases for columns time (in MJD), bandpass Name, fivesigmaDepth
            (mags)
        photParams: instance of `sims.photUtils.PhotParams`, optional,
            defaults to None
            Describes the observing conditions and telescope for the
            observation. The default value of None instantiates this for the
            LSST site and telescope.

        Returns
        -------
        astropy.Table containing the following columns at the minimum
        ['time', 'flux', 'fluxerr', 'band']. If a set of indexes are provided
        for the m5Values
        """
        format = format.lower()

        if format not in ['dataframe']:
            raise NotImplementedError(
                "Unavailable input format {0}".format(format))

        return


if __name__ == '__main__':
    config = bcm.BaseCatalogConfig()
    config.load(os.path.join(getPackageDir("sims_catUtils"), "config",
                             "db.py"))
    username = DbAuth.username(config.host, config.port)
    password = DbAuth.password(config.host, config.port)
    DBConnection = pymssql.connect(user=username,
                                   password=password,
                                   database=config.database,
                                   port=config.port)
    db = DBConnection.cursor()
    reflc = RefLightCurves(idSequence=(6001163623700, 6000324908000),
                           tableName='TwinkSN',
                           dbConnection=DBConnection,
                           dbCursor=db)
    print(reflc)
    print(reflc.columns)
    print(reflc.dbConnection)
    print(reflc.idSequence)
    print('query_stmt', reflc.query())
    print('dataframe', reflc.get_params())
    sn = reflc.astro_object(idValue=6001163623700, obsMetaDataList=None)
    print(sn)
    print(sn.SNstate)
