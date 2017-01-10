"""
Unit tests for Monitor class.
"""
from __future__ import absolute_import
import os
import subprocess
import unittest
import desc.monitor

class MonitorTestCase(unittest.TestCase):
    def setUp(self):
        self.test_db = 'test_cache.db'
        subprocess.call('touch %s' % self.test_db, shell=True)
        self.dbConn = desc.monitor.dbInterface(database=self.test_db,
                                               host=None, port=None,
                                               driver='sqlite')
        self.monitor = desc.monitor.Monitor(self.dbConn)
        self.photfile = 'light_curve.txt'
        subprocess.call('touch %s' % self.photfile, shell=True)
    def tearDown(self):
        os.remove(self.photfile)
        del self.monitor
    def test_db_load(self):
        self.assertEqual(self.monitor.num_visits, 1)

if __name__ == '__main__':
    unittest.main()
