from __future__ import absolute_import
import numpy as np
from lsst.sims.catalogs.generation.db import DBObject

class dbInterface(object):

    def __init__(self, database, host, port, driver):

        self._dbo = DBObject(database=database, host=host, port=port,
                             driver=driver)

    def forcedSourceFromId(self, objectId):

        dtype = np.dtype([('object_id', np.int),
                          ('ccd_visit_id', np.int),
                          ('psf_flux', np.float),
                          ('psf_flux_err', np.float),
                          ('flags', np.int)])
        query = """select objectId, ccdVisitId, psFlux, psFlux_Sigma, flags
                   from ForcedSource
                   where objectId = %i""" % objectId
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def visitFromCcdVisitId(self, visitId):

        dtype = np.dtype([('visit_id', np.int),
                          ('filter', str, 300),
                          ('obs_start', str, 300)])

        query = """select visitId, filterName, obsStart
                   from CcdVisit
                   where ccdVisitId = %i""" % visitId
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results

    def objectFromId(self, objectId):

        dtype = np.dtype([('object_id', np.int),
                          ('parent_object_id', np.int),
                          ('ra', np.float),
                          ('dec', np.float)])
        query = """select objectId, parentObjectId, psRa, psDecl
                   from Object
                   where objectId = %i""" % objectId
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
                   and (psDecl > %f) and (psDecl < %f)""" % (ra - tol,
                                                             ra + tol,
                                                             dec - tol,
                                                             dec + tol)
        results = self._dbo.execute_arbitrary(query, dtype=dtype)
        return results
