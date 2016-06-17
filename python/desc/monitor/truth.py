"""
A class to download the truth parameters of objects and build light curves
from them.
"""

from __future__ import absolute_import, division, print_function
"""
Module to obtain truth values and truth light curves from the catsim database
"""
import pandas as pd
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.supernovae import SNObject
# import pymssql
# import os
# from lsst.utils import getPackageDir
# from lsst.daf.persistence import DbAuth


class RefLightCurves(object):
    """
    Class to connect to the database tables on fatboy and obtain reference
    light curves

    Parameters
    ----------
    tableName: string, mandatory
        case insensitive string name of table on database to connect to
        for model parameters of astrophysical objects
    idCol: string, optional, defaults to 'snid'
        column name of Index on the table
    columns: tuple of strings, optional, defaults to values for SN
        tuple of strings that completely specify the truth values for
        the astrophysical object
    dbConnection: `pymssql.connection` instance, mandatory
        connection to the database where the relevant tables of catsim
        objects are stored
    dbCursor: `pymssql.connection.cursor` instnce, optional, defaults to None
        cursor to the DataBase Connection. If None, a new cursor is obtained
        from self.dbConnection
    idSequence: sequence of one dimension, optional, defaults to None
        sequence of unique ids in the catsim universe indexing the
        astrophysical objects in the database.

    Examples
    --------
    >>> reflc = RefLightCurves(idSequence=(6001163623700, 6000324908000),
                               tableName='TwinkSN',
                               dbConnection=DBConnection,
                               dbCursor=db) # doctest: +SKIP
    """
    def __init__(self,
                 tableName,
                 idCol='snid',
                 columns=('snid', 'redshift', 'snra', 'sndec', 't0', 'x0',
                          'x1', 'c'),
                 dbConnection=None,
                 dbCursor=None,
                 idSequence=None):

        self.columns = columns
        self.dbConnection = dbConnection
        self._dbCursor = dbCursor
        self._idSequence = idSequence
        self.columns = columns
        self.idCol = idCol
        self.tableName = tableName
        self._idvals = None
        self.objectID = 42

    @property
    def dbCursor(self):
        if self._dbCursor is None:
            self._dbCursor = self.dbConnection.cursor()
        return self._dbCursor

    @staticmethod
    def uniqueIDtoTableId(uniqueID, objTypeID, nshift=10):
        id = uniqueID - objID
        return np.right_shift(id, nshift) 

    @property
    def idSequence(self):
        x = np.asarray(self._idSequence)
        return self.uniqueIDtoTableId(x, objTypeID=42, nshift=10)


    def allIdinTable(self, sqlconstraint='', chunksize=None):
        """
        return a `pd.Series`of all IDs in the table with an optional
        constraint specified as sqlconstraint. If chunkSize is not
        None, but set to an integer, it returns a generator to the
        series returning chunkSize values at a time.

        Parameters
        ----------
        sqlconstraint: string, optional, defaults to ''
            sql constraint specified through a WHERE clause
        chunksize: integer, optional, defaults to None
            if not None, the return value is a generator to
            a series getting chunkSize values at a time

        Returns
        -------
        `pd.Series` of snids or a generator to it returning chunkSize values
        at at time

        Examples
        --------
        >>> ids = reflc.allIdinTable(chunksize=None) # doctest: +SKIP
        >>> ids.astype(int).values.flatten() # doctest: +SKIP
        >>> # `numpy.ndarray` of dtype int, having all SNIDs
        >>> ids = reflc.allIdinTable(chunksize=50) # doctest: +SKIP
        >>> idsnext().astype(int).values.flatten() # doctest: +SKIP
        >>> # `numpy.ndarray` of dtype int, having 50 SNIDs,  repeat
        """
        query = """SELECT {0} FROM {1}""".format(self.idCol, self.tableName)
        query += sqlconstraint

        x = pd.read_sql_query(query, con=self.dbConnection,
                              chunksize=chunksize)
        return x

    def get_numObjects(self, sqlconstraint=''):
        """
        return the number of objects in self.table

        Parameters
        ----------
        sqlconstraint : string, optional, defaults to ''
            sql constraint specified through a WHERE clause

        Returns
        -------
        integer number of objects in table (satisfying legal sqlconstraints if
        specified)

        Examples
        --------
        >>> reflc.get_numObjects() #doctest: +SKIP 
        >>> 776620
        """
        query = """SELECT COUNT(*) FROM {} """.format(self.tableName)
        query += sqlconstraint
        self.dbCursor.execute(query)
        n = self.dbCursor.fetchone()[0]
        return n

    def buildquery(self, idValue=None, columns=None):
        """
        Return the query statement to be used to obtain model
        parameters of the set of astrophysical objects

        Parameters
        ----------
        idValue: integer, optional, defaults to None
            unique ID of the astrophysical object. If None, then the query
            is built for all objects in idSequence. If None, then
            the query is built for all objects in the table
        columns: tuple of strings, optional, defaults to None
            Columns that will be queried, if None, defaults to
            self.columns

        Returns
        -------
        String with a query statement to use to obtain parameters
        """
        if columns is None:
            columns = self.columns

        query = """SELECT """
        query += ", ".join(xx.strip() for xx in columns)
        query += " FROM {} ".format(self.tableName)

        # if query is for a single idvalue, construct the query stmt
        if idValue is not None:
            tableIdValue = self.uniqueIDtoTableId(idValue, objID=42, nshift=10)
            query += "WHERE {0} = {1}".format(self.idCol, tableIdValue)
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

        Parameters
        ----------
        idValue: int, optional, defaults to None
            indexes of the astrophysical objects whose parameters are to be
            obtained. If idValue is None, then all astrophysical objects with
            ids in `self.idSequence` are used. If `self.idSequence` is None,
            then all astrophysical objects are used.
        """
        df = pd.read_sql_query(self.buildquery(idValue=idValue),
                               self.dbConnection,
                               index_col=self.idCol,
                               coerce_float=False)
        return df

    def astro_object(self, idValue, mjdOffset=59580.):
        """
        instance of the catsim representation of the astrophysical object.

        Parameters
        ----------
        idValue: int, mandatory
            index of the astro_object
        mjdOffset: float, optional, defaults to 59580.
            offset in time parameters in the database for the transient objects

        Returns
        -------
        Instance of astro_object with parameters from the database

        Examples
        --------
        >>> sn = reflc.astro_object(idValue=6001163623700)

        """
        df = self.get_params(idValue)
        sn = SNObject(ra=df.snra.values[0], dec=df.sndec.values[0])
        paramDict = dict()
        for param in ['t0', 'x0', 'x1', 'c']:
            paramDict[param] = df[param].values
        paramDict['t0'] += mjdOffset
        paramDict['z'] = df.redshift.values[0]
        sn.set(**paramDict)
        return sn

    def lightCurve(self, idValue, observations, bandpassDict,
                   format='dataframe',
                   keys=['expMJD', 'filter', 'fiveSigmaDepth'],
                   photParams=None):
        """
        return the light curve of the object for observations as a
        `pandas.Dataframe`.
        Parameters
        ----------
        idValue: integer, mandatory
            index of astrophysical object
        observations: mandatory, allowed formats decided by format variable
            observations corresponding to which the light curve is obtained
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
        astropy.Table containing the following columns at the minimum
        ['time', 'flux', 'fluxerr', 'band']. If a set of indexes are provided
        for the m5Values (eg. as indexes in a dataFrame)
        """
        format = format.lower()

        if format not in ['dataframe']:
            raise NotImplementedError(
                "Unavailable input format {0}".format(format))

        sn = self.astro_object(idValue)
        timeMin = sn.mintime()
        timeMax = sn.maxtime()

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
                                     bandpassobject=bandpassDict[band])
            m5val = df.ix[ind, 'fiveSigmaDepth']
            fluxerr = sn.catsimBandFluxError(time=time,
                                             bandpassobject=bandpassDict[band],
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
        return output
