from __future__ import absolute_import
import numpy as np
from lsst.sims.catalogs.db import DBObject

__all__ = ['dbInterface']

class dbInterface(object):

    def __init__(self, database, host, port, driver,
                 project='Twinkles Run1.1'):

        self._dbo = DBObject(database=database, host=host, port=port,
                             driver=driver)
        self.project = project

    def forcedSourceFromId(self, objectId):

        dtype = np.dtype([('object_id', np.int),
                          ('ccd_visit_id', np.int),
                          ('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int)])
        query = """select objectId, ccdVisitId, psFlux, psFlux_Sigma, flags
                   from ForcedSource
                   where objectId = %i
                   and project = %s""" % (objectId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def visitFromCcdVisitId(self, visitId):

        dtype = np.dtype([('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300)])

        query = """select visitId, filterName, obsStart
                   from CcdVisit
                   where ccdVisitId = %i
                   and project = %s""" % (visitId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def all_fs_visits_from_id(self, objectId):

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
                       and ForcedSource.project = %s""" % (objectId,
                                                           self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def objectFromId(self, objectId):

        dtype = np.dtype([('object_id', np.int),
                          ('parent_object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, parentObjectId, psRa, psDecl
                   from Object
                   where objectId = %i
                   and project = %s""" % (objectId, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def objectFromRaDec(self, ra, dec, tol):

        dtype = np.dtype([('object_id', np.int),
                          ('parent_object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, parentObjectId, psRa, psDecl
                   from Object
                   where (psRa > %f) and (psRa < %f)
                   and (psDecl > %f) and (psDecl < %f)
                   and project = %s""" % (ra - tol, ra + tol,
                                          dec - tol, dec + tol, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def get_all_visit_info(self):

        dtype = np.dtype([('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300)])

        query = """select visitId, filterName, obsStart
                   from CcdVisit where project = %s""" % self.project
        results = self._dbo.execute_arbitrary(query, dtype=dtype)

        return results

    def get_all_objects_in_visit(self, ccd_visit_id):

        dtype = np.dtype([('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int)])
        query = """select psFlux, psFlux_Sigma, flags
                   from ForcedSource
                   where ccd_visit_id = %i
                   and project = %s""" % (ccd_visit_id, self.project)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results
