"""
Copyright (C) 2016, LSST Dark Energy Science Collaboration (DESC)
All rights reserved.

This software may be modified and distributed under the terms
of the modified BSD license.  See the LICENSE file for details.
"""
# ==============================================================================

from __future__ import print_function

import os
# import pandas as pd
import numpy as np
import sncosmo
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from lsst.utils import getPackageDir

# ==============================================================================

class Monitor(object):
    '''
    Simple class for extracting DM forced photometry light curves.
    '''
    def __init__(self):
        return

    def read(self, photfile):
        print("Cannot read data from", photfile, "yet.")
        return

    def run(self, algorithm=None):
        if algorithm is None:
            print("No algorithms coded yet.")
        return

# =============================================================================

class LightCurve(object):
    '''
    Fetch ForcedSource data with Butler and package for use with accompanying
    visualization routines.
    '''
    def __init__(self, fp_table_dir=None, mjd_file=None):

        self.fp_table_dir = fp_table_dir
        self.visit_map = {}
        self.bandpasses = []
        self.lightcurve = None
        self.visit_mjd = {}

        # Associate mjd with visits (just a hack for now, need to think about this more):
        mjd_list = np.genfromtxt(mjd_file, delimiter=',')
        for visit_date in mjd_list:
            self.visit_mjd[str(int(visit_date[0]))] = visit_date[1]

        for visit_dir in os.listdir(str(fp_table_dir+'/0/')):
            visit_band = visit_dir[-1]
            visit_num = visit_dir[1:-3]
            if visit_band not in self.bandpasses:
                self.bandpasses.append(visit_band)
                self.visit_map[visit_band] = [visit_num]
                # Register required lsst bandpass in sncosmo registry:
                bandpass_file = os.path.join(str(getPackageDir('throughputs') +
                                                 '/baseline/total_' +
                                                 visit_band + '.dat'))
                bandpass_info = np.genfromtxt(bandpass_file,
                                              names=['wavelen', 'transmission'])
                band = sncosmo.Bandpass(bandpass_info['wavelen'],
                                        bandpass_info['transmission'],
                                        name=str('lsst' + visit_band),
                                        wave_unit=u.nm)
                try:
                    sncosmo.registry.register(band)
                except Exception:
                    continue
            else:
                self.visit_map[visit_band].append(visit_num)

    def build_lightcurve(self, objid):
        """
        Assemble a light curve data table from available files.
        """
        lightcurve = {}
        lightcurve['bandpass'] = []
        lightcurve['mjd'] = []
        lightcurve['ra'] = []
        lightcurve['dec'] = []
        lightcurve['flux'] = []
        lightcurve['flux_error'] = []
        lightcurve['zp'] = []
        lightcurve['zpsys'] = []

        for bandpass in self.bandpasses:
            for visit in self.visit_map[bandpass]:
                # NB. Hard-coded filename conventions:
                hdulist = fits.open(str(self.fp_table_dir + '/0/v' + str(visit) +
                                        '-f'+bandpass+'/R22/S11.fits'))
                obj_data = hdulist[1].data[np.where(hdulist[1].data['objectId'] == objid)]
                if len(obj_data) > 0:
                    lightcurve['bandpass'].append(str('lsst' + bandpass))
                    lightcurve['mjd'].append(self.visit_mjd[str(visit)])
                    lightcurve['ra'].append(obj_data['coord_ra'][0])
                    lightcurve['dec'].append(obj_data['coord_dec'][0])
                    lightcurve['flux'].append(obj_data['base_PsfFlux_flux'][0])
                    lightcurve['flux_error'].append(obj_data['base_PsfFlux_fluxSigma'][0])
                    lightcurve['zp'].append(25.0)
                    lightcurve['zpsys'].append('ab')
        self.lightcurve = Table(data=lightcurve)


    def visualize_lightcurve(self):
        """
        Make a simple light curve plot.
        """
        if self.lightcurve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        fig = sncosmo.plot_lc(self.lightcurve)

        return fig

# ==============================================================================
