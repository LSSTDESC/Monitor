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
        self.monitor = desc.monitor.Monitor()
        self.photfile = 'light_curve.txt'
        subprocess.call('touch %s' % self.photfile, shell=True)
    def tearDown(self):
        os.remove(self.photfile)
        del self.monitor
    def test_read(self):
        self.monitor.read(self.photfile)
    def test_run(self):
        self.monitor.read(self.photfile)
        self.monitor.run()

if __name__ == '__main__':
    unittest.main()
