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
import pandas as pd
import numpy as np
import sncosmo
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from lsst.utils import getPackageDir
plt.style.use('ggplot')

__all__ = ['Monitor', 'LightCurve', 'SeeingCurve']
# ==============================================================================

class Monitor(object):
    '''
    Simple class for extracting DM forced photometry light curves.
    '''
    def __init__(self, dbConn, truth_dbConn=None):
        self.dbConn = dbConn
        self.truth_dbConn = truth_dbConn
        self.return_lightcurve = {}
        self.num_visits = self.dbConn.get_number_of_visits()
        self.best_seeing = None
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

        if self.best_seeing is None:
            self.get_best_seeing_visit()

        visit_data = self.dbConn.get_all_visit_info()
        stars = self.get_stars(in_visit=self.best_seeing['ccd_visit_id'],
                               with_sigma_clipping=True)
        visit_flux_err = []

        for star in stars:
            star_data = self.dbConn.all_fs_visits_from_id(star['object_id'])
            self.star_data = star_data
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
        dc.build_lightcurve(depth_curve)

        return dc

    def measure_seeing_curve(self):
        """
        Find the seeing in each visit.
        """

        visit_data = self.dbConn.get_all_visit_info()
        seeing_info = {}
        seeing_info['bandpass'] = [str('lsst' + x) for x in visit_data['filter']]
        timestamp = Time(visit_data['obs_start'], scale='utc')
        seeing_info['mjd'] = timestamp.mjd
        seeing_info['obsHistId'] = visit_data['visit_id']
        seeing_info['seeing'] = visit_data['seeing']

        sc = SeeingCurve(self.dbConn)
        sc.build_seeing_curve(seeing_info)

        return sc

    def get_stars(self, in_visit=None, fainter_than=None, with_sigma_clipping=False):
        """
        Get the stars from the pserv database for a visit.

        Can return only stars fainter than a given magnitude.
        """
        visit_data = self.dbConn.get_all_objects_in_visit(in_visit)

        if fainter_than is not None:
            max_flux = np.power(10, (fainter_than - 22.5)/-2.5)
            visit_data = visit_data[np.where(visit_data['psf_flux'] < max_flux)]

        if with_sigma_clipping is True:
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
        else:
            test_objects = visit_data

        return test_objects

    def get_best_seeing_visit(self):
        """
        Get the i-band visit with the best seeing in arcseconds.
        """

        visit_data = self.dbConn.get_all_visit_info()
        i_visits = visit_data[np.where(visit_data['filter']=='i')]
        best_i_visit = i_visits[np.argmin(i_visits['seeing'])]
        self.best_seeing = best_i_visit




# =============================================================================


class BaseCurve(object):
    """
    A Base Class used to initialize the curve methods: LightCurve and
    SeeingCurve.
    """

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

class LightCurve(BaseCurve):

    def build_lightcurve(self, lc_dict):
        """
        A wrapper around pd.read_dict to build a lightcurve from a dict.
        """

        self.lightcurve = pd.DataFrame.from_dict(lc_dict)

    def build_lightcurve_from_fp_table(self, objid):
        """
        Assemble a light curve data table from available forced photometry
        in a data repo.
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
        self.build_lightcurve(lightcurve)

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

        self.build_lightcurve(lightcurve)

    def visualize_lightcurve(self, using='flux', include_errors=True,
                             use_existing_fig = None):
        """
        Make a simple light curve plot.

        Adapted from sncosmo.plot_lc source code.
        """
        if self.lightcurve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

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
            filt_mjd = self.lightcurve['mjd'][self.lightcurve['bandpass']==filt_name].values
            using_mjd = self.lightcurve[using][self.lightcurve['bandpass']==filt_name].values
            if include_errors is True:
                using_error_mjd = self.lightcurve[str(using+'_error')][self.lightcurve['bandpass']==filt_name].values
                plt.errorbar(filt_mjd, using_mjd, yerr=using_error_mjd,
                         ls='None', marker='.', ms=3, c=color[plot_num-1])
            else:
                plt.scatter(filt_mjd, using_mjd, c='purple', marker='+')
            plt.locator_params(axis='x',nbins=5)
            plt.xlabel('MJD')
            if using == 'flux':
                plt.ylabel(str(using + ' (nmgy)'))
                plt.ylim(-10, np.max(using_mjd)+10)
            elif using == 'mag':
                plt.ylabel(using)
                if use_existing_fig is not None:
                    plt.gca().invert_yaxis()
            plot_num += 1
        plt.tight_layout()

        return fig

# ==============================================================================

class SeeingCurve(BaseCurve):
    """
    An object to store the seeing information from all visits in the survey.
    """

    def build_seeing_curve(self, sc_dict):
        """
        A wrapper around pd.read_dict to build a lightcurve from a dict.
        """

        self.seeing_curve = pd.DataFrame.from_dict(sc_dict)

    def visualize_seeing_curve(self):

        if self.seeing_curve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        n_subplot = len(self.filter_list)*2
        n_col = 2
        n_row = len(self.filter_list)
        fig = plt.figure(figsize = (4. * n_col, 3. * n_row))

        color = ['b', 'g', 'y', 'orange', 'r', 'k']

        plot_num = 1
        for filt in self.filter_list:

            fig.add_subplot(n_row, n_col, plot_num)
            filt_name = str('lsst' + filt)
            plt.title(filt_name)

            filt_seeing = self.seeing_curve['seeing'][self.seeing_curve['bandpass'] == filt_name]
            filt_mjd = self.seeing_curve['mjd'][self.seeing_curve['bandpass'] == filt_name]
            plt.scatter(filt_mjd, filt_seeing,
                        c = color[(plot_num - 1)//2], marker = '+')
            plt.locator_params(axis='x',nbins=5)
            plt.ylabel('Seeing (arcseconds)')
            plt.xlabel('MJD')
            plot_num += 1

            fig.add_subplot(n_row, n_col, plot_num)
            plt.title(filt_name)
            plt.hist(filt_seeing, color=color[(plot_num - 1)//2])
            plt.xlabel('Seeing (arcseconds)')
            plot_num += 1

        fig.tight_layout()

        return fig
