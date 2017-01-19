"""
Unit tests for Monitor class.
"""
from __future__ import absolute_import
import os
import subprocess
import unittest
import sqlite3
import desc.monitor

class MonitorTestCase(unittest.TestCase):
    # def setUp(self):
    #     self.test_db = 'test_cache.db'
    #     self.dbConn = desc.monitor.dbInterface(database=self.test_db,
    #                                            host=None, port=None,
    #                                            driver='sqlite')
    #     self.monitor = desc.monitor.Monitor(self.dbConn)
    #     self.photfile = 'light_curve.txt'
    #     subprocess.call('touch %s' % self.photfile, shell=True)
    # def tearDown(self):
    #     os.remove(self.photfile)
    #     del self.monitor

    def setUp(self):
        """
        Create a test dataset as an sqlite db for exercising
        get_twinkles_vists.
        """
        self._test_file = 'my_test_sqlite.db'
        connection = sqlite3.connect(self._test_file)
        connection.execute("create table CcdVisit (ccdVisitId, project)")
        self._obsHistIDs = (359, 123, 59025, 430201, 430201)
        self._projectID = "'Twinkles Run1.1'"
        for obsHistID in self._obsHistIDs:
            connection.execute("""insert into CcdVisit values
                                  (%i, %s)""" % (obsHistID, self._projectID))
        other_projectID = "'Twinkles Run3'"
        other_obsHistIDs = (340298, 2342, 1550)
        for obsHistID in other_obsHistIDs:
            connection.execute("""insert into CcdVisit values
                                  (%i, %s)""" % (obsHistID, other_projectID))
        connection.commit()
        connection.close()

        self.dbConn = desc.monitor.dbInterface(database=self._test_file,
                                               host=None, port=None,
                                               driver='sqlite')
        self.monitor = desc.monitor.Monitor(self.dbConn)

    def tearDown(self):
        "Clean up test file."
        os.remove(self._test_file)

    def test_db_load(self):
        self.assertEqual(self.monitor.num_visits, 5)

if __name__ == '__main__':
    unittest.main()
