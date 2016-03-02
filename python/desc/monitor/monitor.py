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
    def __init__(self, fp_table_dir=None, bandpasses=None, visit_lists=None, mjd_file=None):

        self.fp_table_dir = fp_table_dir
        self.visit_map = {}
        self.bandpasses = bandpasses
        self.lightcurve = None
        self.visit_mjd = {}

        # Associate mjd with visits (just a hack for now, need to think about this more):
        mjd_list = np.genfromtxt(mjd_file, delimiter=',')
        for visit_date in mjd_list:
            self.visit_mjd[str(int(visit_date[0]))] = visit_date[1]

        for bandpass, visit_list in zip(self.bandpasses, visit_lists):
            # Associate correct visits with bandpass:
            self.visit_map[bandpass] = visit_list
            # Register required lsst bandpass in sncosmo registry:
            bandpass_file = os.path.join(str(getPackageDir('throughputs') +
                                             '/baseline/total_' +
                                             bandpass + '.dat'))
            bandpass_info = np.genfromtxt(bandpass_file,
                                          names=['wavelen', 'transmission'])
            band = sncosmo.Bandpass(bandpass_info['wavelen'],
                                    bandpass_info['transmission'],
                                    name=str('lsst' + bandpass),
                                    wave_unit=u.nm)
            sncosmo.registry.register(band)


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
                hdulist = fits.open(str(self.fp_table_dir + '/0/v' +
                                        str(visit) + '-fr/R22/S11.fits'))
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
        self.lightcurve = lightcurve


    def visualize_lightcurve(self):
        """
        Make a simple light curve plot.
        """
        if self.lightcurve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        lc_table = Table(data=self.lightcurve)
        fig = sncosmo.plot_lc(lc_table)

        return fig


# ==============================================================================
