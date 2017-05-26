from __future__ import absolute_import
import numpy as np
from lsst.sims.catalogs.db import DBObject

__all__ = ['DBInterface', 'TruthDBInterface', 'OpsimDBInterface']


class DBInterface(object):
    """
    Tools designed to access the pserv database.

    Parameters
    ----------
    database : str
        Name of the database.
    host : str
        The location of the host.
    port : int
        Port to use for the connection.
    driver : str
        Type of database. Probably 'mysql'.
    project : str, default="Twinkles Run1.1"
        The specific project for which to return data from the pserv database.
    """

    def __init__(self, database, host, port, driver,
                 project="Twinkles Run1.1"):

        self._dbo = DBObject(database=database, host=host, port=port,
                             driver=driver)
        self.project = project

    def forcedSourceFromId(self, objectId):
        """
        Get the flux information for all possible visits for a specific source.

        Parameters
        ----------
        objectId : int
            The object id of the source

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('ccd_visit_id', np.int),
                          ('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int)])
        query = """select objectId, ccdVisitId, psFlux, psFlux_Sigma, flags
                   from ForcedSource
                   where objectId = %i
                   and project = '%s'""" % (objectId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def visitFromCcdVisitId(self, visitId):
        """
        Get the bandpass and mjd of a visit.

        Parameters
        ----------
        visitId : int
            The ccd visit id of the visit. Not the same as the obsHistId from
            Opsim.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300)])

        query = """select visitId, filterName, obsStart
                   from CcdVisit
                   where ccdVisitId = %i
                   and project = '%s'""" % (visitId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def all_fs_visits_from_id(self, objectId):
        """
        Get the flux information as well as a join on the visit information
        for all possible visits for a specific source.

        Parameters
        ----------
        objectId : int
            The object id of the source

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('ccd_visit_id', np.int),
                          ('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int),
                          ('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300)])
        query = """
                select ForcedSource.objectId, ForcedSource.ccdVisitId,
                       ForcedSource.psFlux, ForcedSource.psFlux_Sigma,
                       ForcedSource.flags, CcdVisit.visitId,
                       CcdVisit.filterName, CcdVisit.obsStart
                from ForcedSource
                inner join CcdVisit on
                       ForcedSource.ccdVisitId = CcdVisit.ccdVisitId
                       and ForcedSource.objectId = %i
                       and ForcedSource.project = '%s'
                       and CcdVisit.project = '%s'""" % (objectId,
                                                         self.project,
                                                         self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def objectFromId(self, objectId):
        """
        Get the location information and parent id for a specific source.

        Parameters
        ----------
        objectId : int
            The object id of the source

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('parent_object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, parentObjectId, psRa, psDecl
                   from Object
                   where objectId = %i
                   and project = '%s'""" % (objectId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def objectFromRaDec(self, ra, dec, tol):
        """
        Get the list of possible objects in a square box of 2*tol degrees
        around the location (ra, dec) in degrees.

        Parameters
        ----------
        ra : float
            RA location in degrees
        dec : float
            DEC location in degrees
        tol : float
            radius of search in degrees

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('parent_object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, parentObjectId, psRa, psDecl
                   from Object
                   where (psRa > %f) and (psRa < %f)
                   and (psDecl > %f) and (psDecl < %f)
                   and project = '%s'""" % (ra - tol, ra + tol,
                                            dec - tol, dec + tol, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def get_all_visit_info(self):
        """
        Get all the visit information for all possible visits in the project.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('ccd_visit_id', np.int),
                          ('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300),
                          ('zero_point', np.float),
                          ('seeing', np.float),
                          ('sky_noise', np.float)])

        query = """select ccdVisitId, visitId, filterName, obsStart, zeroPoint,
                          seeing, skyNoise
                   from CcdVisit where project = '%s'""" % self.project
        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results

    def get_all_objects_in_visit(self, ccd_visit_id):
        """
        Get the flux and location information for all sources in a visit.

        Parameters
        ----------
        ccd_visit_id : int
            The ccd visit id of the visit. Different from obsHistId from Opsim.
            Can get the ccd visit id for visits using get_all_visit_info.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float),
                          ('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int)])
        query = """select ForcedSource.objectId,
                          Object.psRa, Object.psDecl,
                          ForcedSource.psFlux,
                          ForcedSource.psFlux_Sigma,
                          ForcedSource.flags
                   from ForcedSource
                   inner join Object on
                          ForcedSource.objectId = Object.objectId
                          and ForcedSource.ccdVisitId = %i
                          and ForcedSource.project = '%s'
                          and Object.project = '%s'""" % (ccd_visit_id,
                                                          self.project,
                                                          self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results

    def get_number_of_visits(self):
        """
        Get the count on the ccdVisitId's in the project.

        Returns
        -------
        results : int
            Number of visits in the project.
        """
        dtype = np.dtype([('count', np.int)])

        query = """select count(ccdVisitId)
                   from CcdVisit
                   where project = '%s'""" % (self.project)

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results['count'][0]

    def get_object_locations(self):
        """
        Get the locations for all objects in a project.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query to the pserv database.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, psRa, psDecl
                   from Object
                   where project = '%s'""" % (self.project)

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results


class TruthDBInterface(object):
    """
    Tools designed to access the simulation "truth" database.

    Parameters
    ----------
    database : str
        Name of the database.
    host : str or None, default=None
        The location of the host.
    port : int or None, default=None
        Port to use for the connection.
    driver : str, default='sqlite'
        Type of database. Probably 'sqlite'.
    """

    def __init__(self, database, host=None, port=None, driver='sqlite'):

        self._dbo = DBObject(database=database, host=host, port=port,
                             driver=driver)

    def get_stars_by_visit(self, visit_num):
        """
        Get location and flux for all the stars for a visit.

        Parameters
        ----------
        visit_num : int
            The obsHistId of the visit.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query.
        """

        dtype = np.dtype([('object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float),
                          ('true_flux', np.float),
                          ('true_flux_error', np.float)])
        query = """SELECT unique_id, ra, dec, true_flux, true_flux_error
                   FROM stars
                   WHERE obsHistId = %i""" % (visit_num)

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results

    def get_all_star_ra_dec(self):
        """
        Get the location for all the stars in the database.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query.
        """

        dtype = np.dtype([('object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """SELECT unique_id, ra, dec
                   FROM stars"""

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results

    def get_all_info_by_object_id(self, object_id_tuple):
        """
        Get location and flux for a tuple of object_ids.

        Parameters
        ----------
        object_id_tuple : tuple of ints
            The object_ids of the sources.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query.
        """
        dtype = np.dtype([('object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float),
                          ('filter', str, 300),
                          ('true_flux', np.float),
                          ('true_flux_error', np.float),
                          ('obsHistId', np.int)])

        query = str("""SELECT *
                   FROM stars
                   WHERE """ + ' or '.join(('unique_id = ' +
                                            str(n) for n in object_id_tuple)))

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results


class OpsimDBInterface(object):
    """
    Tools designed to access an LSST Opsim Database.

    Parameters
    ----------
    database : str
        Name of the database.
    host : str or None, default=None
        The location of the host.
    port : int or None, default=None
        Port to use for the connection.
    driver : str, default='sqlite'
        Type of database. Probably 'sqlite'.
    """
    def __init__(self, database, host=None, port=None, driver='sqlite'):

        self._dbo = DBObject(database=database, host=host, port=port,
                             driver=driver)

    def get_summary_depth_info_for_visits(self, fieldID):
        """
        Get visit observing statistics information from Opsim database for all
        visits to a specified field.

        Parameters
        ----------
        fieldID : int
            The id of the observation field.

        Returns
        -------
        results : numpy recarray
            Numpy recarray with the results of the query.
        """

        dtype = np.dtype([('obsHistID', np.int),
                          ('FWHMeff', np.float),
                          ('rawSeeing', np.float),
                          ('fiveSigmaDepth', np.float),
                          ('airmass', np.float),
                          ('filter', str, 300)])

        query = str("""SELECT obsHistID, FWHMeff, rawSeeing, fiveSigmaDepth,
                              airmass, filter
                       FROM Summary
                       WHERE fieldID = %s""" % fieldID)

        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results
