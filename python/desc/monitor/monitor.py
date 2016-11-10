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
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from lsst.utils import getPackageDir
plt.style.use('ggplot')

__all__ = ['Monitor', 'LightCurve']
# ==============================================================================

class Monitor(object):
    '''
    Simple class for extracting DM forced photometry light curves.
    '''
    def __init__(self, dbConn):
        self.dbConn = dbConn
        self.return_lightcurve = {}
        self.num_visits = self.dbConn.get_number_of_visits()
        return

    def get_lightcurves(self, lc_list):
        """
        Get the database information for a list of lightcurves and store it.
        """

        for lc_id in lc_list:
            if lc_id not in self.return_lightcurve.keys():
                self.return_lightcurve[lc_id] = LightCurve(self.dbConn)
                self.return_lightcurve[lc_id].build_lightcurve_from_db(lc_id)

    def measure_depth_curve(self):
        """
        Find the 5-sigma limiting depth in each visit.
        """

        visit_data = self.dbConn.get_all_visit_info()
        stars = self.get_stars()
        visit_flux_err = []
        for star in stars:
            star_data = self.dbConn.all_fs_visits_from_id(star['object_id'])
            visit_flux_err.append(star_data['psf_flux_err'])
        visit_flux_err = np.array(visit_flux_err).T
        visit_depths = np.median(visit_flux_err, axis=1)
        visit_depths = 22.5 - 2.5*np.log10(5*visit_depths)

        depth_curve = {}
        depth_curve['bandpass'] = [str('lsst' + x) for x in visit_data['filter']]
        timestamp = Time(visit_data['obs_start'], scale='utc')
        depth_curve['mjd'] = timestamp.mjd
        depth_curve['obsHistId'] = visit_data['visit_id']
        depth_curve['mag'] = visit_depths
        depth_curve['mag_error'] = np.zeros(self.num_visits)

        dc = LightCurve(self.dbConn)
        dc.lightcurve = Table(data=depth_curve)

        return dc

    def get_stars(self, fainter_than=None):
        """
        Get the stars from the pserv database for a visit.

        Can return only stars fainter than a given magnitude.
        """
        best_visit = self.get_best_seeing_visit()
        visit_data = self.dbConn.get_all_objects_in_visit(best_visit['ccd_visit_id'])

        if fainter_than is not None:
            max_flux = np.power(10, (fainter_than - 22.5)/-2.5)
            visit_data = visit_data[np.where(visit_data['psf_flux'] < max_flux)]

        median_val = np.median(visit_data['psf_flux'])
        std_val = np.std(visit_data['psf_flux'])
        max_cut = median_val + (3. * std_val)
        min_cut = median_val - (3. * std_val)
        keep_objects = visit_data[np.where((visit_data['psf_flux'] < max_cut)
                                         & (visit_data['psf_flux'] > min_cut))]
        np.random.seed(42)
        test_objects = np.random.choice(keep_objects, size=100, replace=False)
        ### Prune those that do not have values in all visits
        full_visits = []
        for star_num in range(len(test_objects)):
            star_visits = self.dbConn.forcedSourceFromId(test_objects[star_num]['object_id'])
            if len(star_visits) == self.num_visits:
                full_visits.append(star_num)
        test_objects = test_objects[full_visits]

        return test_objects

    def get_best_seeing_visit(self):
        """
        Get the i-band visit with the best seeing in arcseconds.
        """

        visit_data = self.dbConn.get_all_visit_info()
        i_visits = visit_data[np.where(visit_data['filter']=='i')]
        best_i_visit = i_visits[np.argmin(i_visits['seeing'])]

        return best_i_visit




# =============================================================================


class LightCurve(object):
    '''
    Fetch ForcedSource data with Butler and package for use with accompanying
    visualization routines.
    '''
    def __init__(self, dbConn, fp_table_dir=None, mjd_file=None,
                 filter_list=['u', 'g', 'r', 'i', 'z', 'y']):

        self.dbConn_lc = dbConn
        self.filter_list = filter_list

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

    def build_lightcurve_from_db(self, objid=None, ra_dec=None,
                                 tol=0.005):
        """
        Assemble a light curve data table from available files.
        """

        if (objid is None) and (ra_dec is None):
            raise ValueError('Must specify either objid or ra,dec location.')
        elif (objid is not None) and (ra_dec is not None):
            raise ValueError(str('Please specify only one of objid or ' +
                                 'ra,dec location.'))

        if objid is not None:
            obj_info = self.dbConn_lc.objectFromId(objid)
            fs_info = self.dbConn_lc.all_fs_visits_from_id(objid)

        if ra_dec is not None:
            obj_info = self.dbConn_lc.objectFromRaDec(ra_dec[0], ra_dec[1], tol)
            if len(obj_info) == 0:
                raise ValueError('No objects within specified ra,dec range.')
            elif len(obj_info) > 1:
                print(obj_info)

        num_results = len(fs_info)
        lightcurve = {}
        lightcurve['bandpass'] = [str('lsst' + x) for x in fs_info['filter']]
        timestamp = Time(fs_info['obs_start'], scale='utc')
        lightcurve['obsHistId'] = fs_info['visit_id']
        lightcurve['mjd'] = timestamp.mjd
        lightcurve['ra'] = [obj_info['ra'][0]]*num_results
        lightcurve['dec'] = [obj_info['dec'][0]]*num_results
        lightcurve['flux'] = fs_info['psf_flux']
        lightcurve['flux_error'] = fs_info['psf_flux_err']
        lightcurve['mag'] = 22.5 - 2.5*np.log10(fs_info['psf_flux'])
        lightcurve['mag_error'] = 2.5*np.log10(1 + (fs_info['psf_flux_err']/fs_info['psf_flux']))
        lightcurve['zp'] = [25.0]*num_results #TEMP
        lightcurve['zpsys'] = ['ab']*num_results #TEMP

        self.lightcurve = Table(data=lightcurve)

    def visualize_lightcurve(self, using='flux', include_errors=True,
                             use_existing_fig = None):
        """
        Make a simple light curve plot.

        Adapted from sncosmo.plot_lc source code.
        """
        if self.lightcurve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        #if using == 'flux':
            #fig = plt.figure()#sncosmo.plot_lc(self.lightcurve)
        n_subplot = len(self.filter_list)
        n_col = 2
        n_row = (n_subplot - 1) // n_col + 1
        if use_existing_fig is None:
            fig = plt.figure(figsize = (4. * n_col, 3. * n_row))
        else:
            fig = use_existing_fig

        color = ['b', 'g', 'y', 'orange', 'r', 'k']

        plot_num = 1
        for filt in self.filter_list:
            fig.add_subplot(n_row, n_col, plot_num)
            filt_name = str('lsst' + filt)
            plt.title(filt_name)
            filt_mjd = self.lightcurve['mjd'][np.where(self.lightcurve['bandpass']==filt_name)]
            using_mjd = self.lightcurve[using][np.where(self.lightcurve['bandpass']==filt_name)]
            if include_errors is True:
                using_error_mjd = self.lightcurve[str(using+'_error')][np.where(self.lightcurve['bandpass']==filt_name)]
                plt.errorbar(filt_mjd, using_mjd, yerr=using_error_mjd,
                         ls='None', marker='.', ms=3, c=color[plot_num-1])
            else:
                plt.scatter(filt_mjd, using_mjd, c='purple', marker='+')
            plt.locator_params(axis='x',nbins=5)
            plt.xlabel('mjd')
            if using == 'flux':
                plt.ylabel(str(using + ' (nmgy)'))
                plt.ylim(-10, np.max(using_mjd)+10)
            elif using == 'mag':
                plt.ylabel(using)
                if use_existing_fig is not None:
                    plt.gca().invert_yaxis()
            plot_num += 1
        plt.tight_layout()
        #elif using == 'mag':
        #    fig = plt.figure()


        return fig

# ==============================================================================
