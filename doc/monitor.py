"""
Copyright (C) 2016, LSST Dark Energy Science Collaboration (DESC)
All rights reserved.

This software may be modified and distributed under the terms
of the modified BSD license.  See the LICENSE file for details.
"""
# ==============================================================================

from __future__ import print_function
from __future__ import absolute_import

import os
# import pandas as pd
import numpy as np
import sncosmo
import astropy.units as u
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from lsst.utils import getPackageDir
from .dbConnection import dbInterface as dbi

__all__ = ['Monitor', 'LightCurve']
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
    def __init__(self, fp_table_dir=None, mjd_file=None,
                 filter_list=['u', 'g', 'r', 'i', 'z', 'y']):

        for filter_name in filter_list:
            bandpass_file = os.path.join(str(getPackageDir('throughputs') +
                                             '/baseline/total_' +
                                             filter_name + '.dat'))
            bandpass_info = np.genfromtxt(bandpass_file,
                                          names=['wavelen', 'transmission'])
            band = sncosmo.Bandpass(bandpass_info['wavelen'],
                                    bandpass_info['transmission'],
                                    name=str('lsst' + filter_name),
                                    wave_unit=u.nm)
            try:
                sncosmo.registry.register(band)
            except Exception:
                continue

        if fp_table_dir is not None:
            self.fp_table_dir = fp_table_dir
            self.bandpasses = filter_list
            self.visit_map = {x:[] for x in self.bandpasses}
            self.lightcurve = None
            self.visit_mjd = {}

            # Associate mjd with visits (just a hack for now, need to think about
            # this more):
            mjd_list = np.genfromtxt(mjd_file, delimiter=',')
            for visit_date in mjd_list:
                self.visit_mjd[str(int(visit_date[0]))] = visit_date[1]

            for visit_dir in os.listdir(str(fp_table_dir+'/0/')):
                visit_band = visit_dir[-1]
                visit_num = visit_dir[1:-3]
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
                hdulist = fits.open(str(self.fp_table_dir + '/0/v' +
                                        str(visit) +
                                        '-f'+bandpass+'/R22/S11.fits'))
                obj_idx = np.where(hdulist[1].data['objectId'] == objid)
                obj_data = hdulist[1].data[obj_idx]
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

    def build_lightcurve_from_db(self, objid=None, ra_dec=None, tol=0.005,
                                 database='DESC_Twinkles_Level_2',
                                 host='127.0.0.1', port='3307',
                                 driver='mysql'):
        """
        Assemble a light curve data table from available files.
        """

        if (objid is None) and (ra_dec is None):
            raise ValueError('Must specify either objid or ra,dec location.')
        elif (objid is not None) and (ra_dec is not None):
            raise ValueError(str('Please specify only one of objid or ' +
                                 'ra,dec location.'))

        dbConn_lc = dbi(database=database, host=host, port=port, driver=driver)

        if objid is not None:
            obj_info = dbConn_lc.objectFromId(objid)
            fs_results = dbConn_lc.forcedSourceFromId(objid)

        if ra_dec is not None:
            obj_info = dbConn_lc.objectFromRaDec(ra_dec[0], ra_dec[1], tol)
            if len(obj_info) == 0:
                raise ValueError('No objects within specified ra,dec range.')
            elif len(obj_info) > 1:
                print(obj_info)

        lightcurve = {}
        lightcurve['bandpass'] = []
        lightcurve['mjd'] = []
        lightcurve['ra'] = []
        lightcurve['dec'] = []
        lightcurve['flux'] = []
        lightcurve['flux_error'] = []
        lightcurve['zp'] = []
        lightcurve['zpsys'] = []

        for visit in fs_results:
            visit_info = dbConn_lc.visitFromCcdVisitId(visit['ccd_visit_id'])
            lightcurve['bandpass'].append(str('lsst' +
                                              visit_info['filter'][0]))
            timestamp = Time(visit_info['obs_start'][0], scale='utc')
            lightcurve['mjd'].append(timestamp.mjd)
            lightcurve['ra'].append(obj_info['ra'][0])
            lightcurve['dec'].append(obj_info['dec'][0])
            lightcurve['flux'].append(visit['psf_flux'])
            lightcurve['flux_error'].append(visit['psf_flux_err'])
            lightcurve['zp'].append(25.0) ### TEMP
            lightcurve['zpsys'].append('ab') ### TEMP

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
