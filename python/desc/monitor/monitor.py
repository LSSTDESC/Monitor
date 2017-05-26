"""
Copyright (C) 2016, LSST Dark Energy Science Collaboration (DESC)
All rights reserved.

This software may be modified and distributed under the terms
of the modified BSD license.  See the LICENSE file for details.
"""
# =============================================================================

from __future__ import print_function
from __future__ import absolute_import, division

import os
import pandas as pd
import numpy as np
import sncosmo
import astropy.units as u
import matplotlib.pyplot as plt
import warnings
from astropy.time import Time
from astropy.io import fits
from lsst.utils import getPackageDir
from astroML.crossmatch import crossmatch_angular
from scipy.stats import sigmaclip
plt.style.use('ggplot')

__all__ = ['Monitor', 'LightCurve', 'SeeingCurve']
# =============================================================================


class Monitor(object):

    '''
    Simple class for extracting DM forced photometry light curves.

    ...

    Parameters
    ----------
    dbConn : dbInterface instance
        This is a connection to a science database with the LSST DM processed
        data products, e.g. a NERSC hosted pserv database.
        For an example see the 'depth_curve_example' notebook
        in the examples folder.
    truth_dbConn : truthDBInterface instance, optional
        This is a connection to a database containing simulation inputs. This
        will be necessary for methods comparing DM processed results to
        simulation inputs such as the match_catalogs method.

    Attributes
    ----------

    return_lightcurve : dictionary
        A dictionary of LightCurve objects where the keys are the ids of the
        objects in the dbConn database.
    num_visits : int
        The total number of visits for which there are objects in the dbConn
        database.
    best_seeing : int or None (default None)
        The visit id of the best seeing i-band visit in the dbConn database.
    matched_ids : pandas dataframe or None (default None)
        A pandas dataframe containing the dbConn database id numbers for
        objects that are matched to the truth database with match_catalogs
        along with the respective truth database object ids.
    flux_stats : pandas dataframe
        Set when using the calc_flux_residual method. Contains information
        on differences between object fluxes in the simulated input "truth"
        catalog and the DM processed outputs in the dbConn database.
    '''

    def __init__(self, dbConn, truth_dbConn=None):

        self.dbConn = dbConn
        self.truth_dbConn = truth_dbConn
        self.return_lightcurve = {}
        self.num_visits = self.dbConn.get_number_of_visits()

        self.best_seeing = None
        self.matched_ids = None
        self.flux_stats = None
        return

    def get_lightcurves(self, lc_list):

        """
        Get the database information for a list of lightcurves and store it.

        Parameters
        ----------
        lc_list : list
            List of object ids for which to get lightcurves from the results
            database.
        """

        for lc_id in lc_list:
            if lc_id not in self.return_lightcurve.keys():
                self.return_lightcurve[lc_id] = LightCurve(self.dbConn)
                self.return_lightcurve[lc_id].build_lightcurve_from_db(lc_id)

    def measure_depth_curve(self, using='DM_modified'):
        """
        Find the 5-sigma limiting depth in each visit and return as lightcurve.

        Parameters
        ----------
        method : str, ('DM_modified', 'DM', 'stars'), default='DM_modified':
            Sets the method for calculating the 5 sigma depth of a visit. The
            'DM' and 'DM_modified' methods use the seeing and sky_noise values
            from the visit table in the results database. See the
            simple_error_model notebook in the examples folder for more
            information. The 'stars' method uses the flux errors from the
            stars in the visits of a Twinkles run to calculate image depth.

        Returns
        -------
        dc : lightcurve object
            The 5-sigma limiting depth at each visit stored in a lightcurve
            object.
        """

        if self.best_seeing is None:
            self.get_best_seeing_visit()

        visit_data = self.dbConn.get_all_visit_info()

        if using == 'stars':

            stars = self.get_stars(in_visit=self.best_seeing['ccd_visit_id'],
                                   with_sigma_clipping=True)
            visit_flux_err = []

            for star in stars:
                star_db = self.dbConn.all_fs_visits_from_id(star['object_id'])
                visit_flux_err.append(star_db['psf_flux_err'])

            visit_flux_err = np.array(visit_flux_err).T
            visit_depths = np.median(visit_flux_err, axis=1)
            visit_depths = 22.5 - 2.5*np.log10(5*visit_depths)

        elif using == 'DM':

            psf_area = (np.pi*(((visit_data['seeing']/2)*(1/.2))**2.))
            sky_noise_5_sigma = (5*1e9*visit_data['sky_noise'])
            visit_depths = 22.5 - 2.5*np.log10((sky_noise_5_sigma * psf_area) /
                                               visit_data['zero_point'])

        elif using == 'DM_modified':

            psf_area = (np.pi*(((visit_data['seeing']/2)*(1/.2))**2.))
            sky_noise = visit_data['sky_noise']*np.sqrt(psf_area)
            f5_signal = ((1e9*(25 + 25*((1+(4*(sky_noise**2)/25))**.5))/2) /
                         visit_data['zero_point'])
            visit_depths = 22.5 - 2.5*np.log10(f5_signal)

        else:
            raise IOError("Don't understand using == %s. Please specify"
                          " 'DM' or 'DM_modified' or 'stars'" % using)

        depth_curve = {}
        depth_curve['bandpass'] = [str('lsst' + x) for x in
                                   visit_data['filter']]
        timestamp = Time(visit_data['obs_start'], scale='utc')
        depth_curve['mjd'] = timestamp.mjd
        depth_curve['obsHistId'] = visit_data['visit_id']
        depth_curve['mag'] = visit_depths
        depth_curve['mag_error'] = np.zeros(self.num_visits)

        dc = LightCurve(self.dbConn)
        dc._build_lightcurve(depth_curve)

        return dc

    def measure_seeing_curve(self):

        """
        Find the seeing in each visit.

        Returns
        -------
        sc : SeeingCurve object
            Stores the seeing for all the visits.
        """

        visit_data = self.dbConn.get_all_visit_info()
        seeing_info = {}
        seeing_info['bandpass'] = [str('lsst' + x) for x in
                                   visit_data['filter']]
        timestamp = Time(visit_data['obs_start'], scale='utc')
        seeing_info['mjd'] = timestamp.mjd
        seeing_info['obsHistId'] = visit_data['visit_id']
        seeing_info['seeing'] = visit_data['seeing']

        sc = SeeingCurve(self.dbConn)
        sc._build_seeing_curve(seeing_info)

        return sc

    def get_stars(self, in_visit, fainter_than=None,
                  with_sigma_clipping=False):
        """
        Get the stars from the pserv database for a visit.

        Can return only stars fainter than a given magnitude.

        Parameters
        ----------
        in_visit : int
            Specifies the visit number of the observations.
        fainter_than : float or None
            Sets the upper limit on magnitude of the stars that should be
            returned.
        with_sigma_clipping : Boolean, default=False
            If set to True then only stars with magnitudes within 3 standard
            deviations of the median value will be returned.

        Returns
        -------
        test_objects : numpy recarray
            Numpy recarray with the stars from the desired visit.
        """
        visit_data = self.dbConn.get_all_objects_in_visit(in_visit)

        if fainter_than is not None:
            max_flux = np.power(10, (fainter_than - 22.5)/-2.5)
            visit_data = visit_data[np.where(visit_data['psf_flux'] <
                                             max_flux)]

        if with_sigma_clipping is True:
            median_val = np.median(visit_data['psf_flux'])
            std_val = np.std(visit_data['psf_flux'])
            max_cut = median_val + (3. * std_val)
            min_cut = median_val - (3. * std_val)
            keep_objects = visit_data[np.where((visit_data['psf_flux'] <
                                                max_cut) &
                                               (visit_data['psf_flux'] >
                                                min_cut))]
            np.random.seed(42)
            test_objects = np.random.choice(keep_objects, size=100,
                                            replace=False)
            # Prune those that do not have values in all visits
            full_visits = []
            for star_num in range(len(test_objects)):
                star_visits = self.dbConn.forcedSourceFromId(
                                           test_objects[star_num]['object_id'])
                if len(star_visits) == self.num_visits:
                    full_visits.append(star_num)
            test_objects = test_objects[full_visits]
        else:
            test_objects = visit_data

        return test_objects

    def get_best_seeing_visit(self):
        """
        Get the i-band visit with the best seeing in arcseconds.

        Attributes
        ----------
        best_seeing : int or None (default None)
            The visit id of the best seeing i-band visit in
            the dbConn database.
        """

        visit_data = self.dbConn.get_all_visit_info()
        i_visits = visit_data[np.where(visit_data['filter'] == 'i')]
        best_i_visit = i_visits[np.argmin(i_visits['seeing'])]
        self.best_seeing = best_i_visit

    def match_catalogs(self, within_radius=1./3600., return_distance=False):
        """
        Match the ids from the truth catalog and the data catalog by
        ra,dec positions.

        Parameters
        ----------
        within_radius : float, default=1./3600.
            The search radius (in degrees) for matching objects in the two
            catalogs. The default setting is one arcsec.
        return_distance : boolean, default=False
            If true this will return the distances of each matched object in
            the data catalog from the object in the truth catalog.

        Returns
        -------
        matched_distance : numpy array, optional
            Returned when return_distance is True. Records distances in
            degrees for the matched data catalog object to the respective
            truth catalog object.

        Attributes
        ----------
        matched_ids : pandas dataframe or None (default None)
            A pandas dataframe containing the dbConn database id numbers for
            objects that are matched to the truth database with match_catalogs
            along with the respective truth database object ids.
        """

        print('Querying Truth Catalog')
        truth_data = self.truth_dbConn.get_all_star_ra_dec()
        print('Querying Pserv Database')
        obs_data = self.dbConn.get_object_locations()

        true_pos = truth_data[['ra', 'dec']]
        obs_pos = obs_data[['ra', 'dec']]

        true_ra_dec = np.empty((len(true_pos), 2), dtype=np.float64)
        obs_ra_dec = np.empty((len(obs_pos), 2), dtype=np.float64)

        true_ra_dec[:, 0] = true_pos['ra']
        true_ra_dec[:, 1] = true_pos['dec']
        obs_ra_dec[:, 0] = obs_pos['ra']
        obs_ra_dec[:, 1] = obs_pos['dec']

        print('Starting angular crossmatch')
        dist, ind = crossmatch_angular(obs_ra_dec, true_ra_dec, within_radius)
        true_obs_match_ids = []

        if return_distance is True:
            matched_distance = []

        for idx in range(len(true_ra_dec)):
            matches = np.where(ind == idx)[0]
            distances = dist[matches]
            match_info = []
            if len(distances) > 0:
                match_info.append(truth_data['object_id'][idx])
                keep = matches[np.argmin(distances)]
                match_info.append(obs_data['object_id'][keep])
                match_info.append(len(distances))
                true_obs_match_ids.append(match_info)
                if return_distance is True:
                    matched_distance.append(np.min(distances))
        print('Done')

        self.matched_ids = pd.DataFrame(true_obs_match_ids,
                                        columns=('truth_object_id',
                                                 'obs_object_id',
                                                 'num_matches'))

        if return_distance is True:
            return np.array(matched_distance)
        else:
            return

    def calc_flux_residuals(self, depth_curve, seeing_curve, for_visits=None):
        """
        Get the truth data and the observed data then subtract the fluxes
        to get the differences.

        Parameters
        ----------
        depth_curve : LightCurve Object
            Depth curve from measure_depth_curve method.
        seeing_curve : SeeingCurve Object
            Seeing curve from measure_seeing_curve method.
        for_visits : list or None, default=None
            List of the visits to use in the comparison. If set to None
            then it uses all available visits.

        Attributes
        ----------
        flux_stats : pandas dataframe
            Contains information on differences between object fluxes in the
            simulated input "truth" catalog and the DM processed outputs in
            the dbConn database.

        """

        dc = depth_curve.lightcurve
        sc = seeing_curve.seeing_curve

        warnings.simplefilter('always', UserWarning)

        if self.matched_ids is None:
            print('Matching Catalogs')
            self.match_catalogs()

        match_dict = dict((k, v) for (k, v) in zip(
                                        self.matched_ids['truth_object_id'],
                                        self.matched_ids['obs_object_id']))

        truth_data = self.truth_dbConn.get_all_info_by_object_id(
                                        self.matched_ids['truth_object_id'])
        truth_data = pd.DataFrame(truth_data)
        id_match = []
        for true_id in truth_data['object_id']:
            id_match.append(match_dict[true_id])
        truth_data['obs_object_id'] = id_match

        print('Gathering Visit Data')
        all_visits = self.dbConn.get_all_visit_info()
        if for_visits is not None:
            ccd_visit_list = []
            visit_list = []
            for visit_num in for_visits:
                ccd_visit = all_visits[np.where(all_visits['visit_id'] ==
                                                visit_num)]['ccd_visit_id']
                if len(ccd_visit) < 1:
                    warnings.warn("You have specified a visit that does" +
                                  " not exist in your project data.  " +
                                  "Skipping over visit %i." % visit_num)
                else:
                    ccd_visit_list.append(ccd_visit[0])
                    visit_list.append(visit_num)
        else:
            ccd_visit_list = all_visits['ccd_visit_id']
            visit_list = all_visits['visit_id']

        print('Querying for object fluxes')
        obs_flux_table = pd.DataFrame(columns=('obs_object_id',
                                               'visit_id',
                                               'filter',
                                               'psf_flux',
                                               'psf_flux_err'))
        obj_queried = 0
        for obs_id in self.matched_ids['obs_object_id']:
            if obj_queried % 100 == 0:
                print('Loaded %i out of %i objects' % (obj_queried,
                                                       len(self.matched_ids)))
            obj_queried += 1
            obj_flux = self.dbConn.all_fs_visits_from_id(obs_id)
            obs_flux_table = obs_flux_table.append(pd.DataFrame(np.array(
                                            [obj_flux['object_id'],
                                             obj_flux['visit_id'],
                                             obj_flux['filter'],
                                             obj_flux['psf_flux'],
                                             obj_flux['psf_flux_err']]).T,
                                            columns=obs_flux_table.columns),
                                            ignore_index=True)

        obs_flux_table['obs_object_id'] = pd.to_numeric(
            obs_flux_table['obs_object_id'])
        obs_flux_table['visit_id'] = pd.to_numeric(
            obs_flux_table['visit_id'])
        obs_flux_table['psf_flux'] = pd.to_numeric(
            obs_flux_table['psf_flux'])
        obs_flux_table['psf_flux_err'] = pd.to_numeric(
            obs_flux_table['psf_flux_err'])

        flux_statistics = []
        for visit_num in visit_list:
            flux_diffs = []
            obs_rows = (obs_flux_table['visit_id'] == visit_num)
            visit_obs = obs_flux_table[['obs_object_id', 'psf_flux']][obs_rows]
            true_rows = (truth_data['obsHistId'] == visit_num)
            visit_true = truth_data[['obs_object_id', 'true_flux']][true_rows]
            visit_data = pd.merge(visit_obs, visit_true, on='obs_object_id')
            flux_diffs = visit_data['psf_flux'] - visit_data['true_flux']
            f_clip, f_low, f_high = sigmaclip(flux_diffs)
            visit_depth = dc['mag'][np.where(dc['obsHistId'] ==
                                             visit_num)[0][0]]
            visit_seeing = sc['seeing'][np.where(sc['obsHistId'] ==
                                                 visit_num)[0][0]]
            visit_bandpass = dc['bandpass'][np.where(dc['obsHistId'] ==
                                                     visit_num)[0][0]]
            flux_statistics.append([visit_num, np.mean(f_clip),
                                    np.mean(f_clip**2.),
                                    visit_depth, visit_seeing,
                                    visit_bandpass, len(flux_diffs)])

        self.flux_stats = pd.DataFrame(flux_statistics,
                                       columns=[('visit_id'), ('mean_resid'),
                                                ('mean_sq_resid'),
                                                ('depth'), ('seeing'),
                                                ('bandpass'), ('num_objects')])

    def _set_stats_plot_limits(self, num_bins=20.):
        """
        Set limits for plots of bias and sigma.

        Parameters
        ----------
        num_bins : float or length 2 list, default=20.
            The number of bins in the depth and seeing respectively. Can use a
            single number to set both bins to the same number.

        Returns
        -------
        p_lims : dict
            Dictionary of maximum and minimum depths and seeing values for
            plotting methods.
        n_bins : length 2 list
            Number of bins in depth and seeing for plotting methods.
        """
        if self.flux_stats is None:
            raise AttributeError("Need to calculate flux_residuals first. " +
                                 "Use self.calc_flux_residuals.")


        if num_bins is not list:
            n_bins = [num_bins, num_bins]

        p_lims = {}
        p_lims['max_d'] = np.max(self.flux_stats['depth'])
        p_lims['min_d'] = np.min(self.flux_stats['depth'])
        p_lims['max_s'] = np.max(self.flux_stats['seeing'])
        p_lims['min_s'] = np.min(self.flux_stats['seeing'])

        return p_lims, n_bins

    def plot_bias_map(self, in_band='r', with_bins=20., use_existing_fig=None,
                      normalize=False):
        """
        Plot the mean residuals of each visit in a 2-d grid as a function of
        depth and seeing.

        Parameters
        ----------
        in_band : str, default = 'r'
            The bandpass for the comparison. ('u', 'g', 'r', 'i', 'z', or 'y')
        with_bins : float or length 2 list, default = 20.
            The number of bins to use to plot the depth and seeing values
            respectively.
        use_existing_fig : matplotlib figure or None, default=None
            Can use an existing matplotlib figure or if set to None will
            return a new figure.
        normalize : Boolean, default=False
            If true, the plot values will be normalized by the standard
            deviation of the flux results in each bin.

        Returns
        -------
        fig : matplotlib figure
            Bias map with depth and seeing as axes.
        """

        p_lims, num_bins = self._set_stats_plot_limits(num_bins=with_bins)
        in_band = str('lsst'+in_band)

        d_grid, s_grid = np.meshgrid(np.linspace(p_lims['min_d']-.01,
                                                 p_lims['max_d']+.01,
                                                 num=num_bins[0]+1),
                                     np.linspace(p_lims['min_s']-.01,
                                                 p_lims['max_s']+.01,
                                                 num=num_bins[1]+1))
        p_lims['delta_d'] = d_grid[0, 1] - d_grid[0, 0]
        p_lims['delta_s'] = s_grid[1, 0] - s_grid[0, 0]
        p_val = np.zeros(np.shape(d_grid))
        for i in range(np.shape(d_grid)[0]):
            for j in range(np.shape(s_grid)[0]):
                depth_vals = self.flux_stats[((self.flux_stats['depth'] >=
                                               d_grid[i, j]) &
                                              (self.flux_stats['depth'] <
                                               (d_grid[i, j] +
                                                p_lims['delta_d'])) &
                                              (self.flux_stats['bandpass'] ==
                                               in_band))]
                grid_vals = depth_vals[((depth_vals['seeing'] >=
                                         s_grid[i, j]) &
                                        (depth_vals['seeing'] <
                                         (s_grid[i, j] + p_lims['delta_s'])))]

                if len(grid_vals) > 0:
                    if normalize is True:
                        p_val[i, j] = np.average(grid_vals['mean_resid'] /
                                                 np.sqrt(
                                                  grid_vals['mean_sq_resid'] -
                                                  grid_vals['mean_resid']**2.),
                                                 weights=grid_vals[
                                                                'num_objects'])
                    else:
                        p_val[i, j] = np.average(grid_vals['mean_resid'],
                                                 weights=grid_vals[
                                                                'num_objects'])

        if use_existing_fig is not None:
            fig = use_existing_fig
        else:
            fig = plt.figure()
        plt.gca().set_axis_bgcolor('k')
        if normalize is True:
            plt.pcolor(d_grid, s_grid, p_val, cmap=plt.cm.coolwarm,
                       vmin=-1., vmax=1.)
        else:
            plt.pcolor(d_grid, s_grid, p_val, cmap=plt.cm.coolwarm,
                       vmin=-0.3, vmax=0.3)
        plt.colorbar()
        plt.xlim(p_lims['min_d'], p_lims['max_d'])
        plt.ylim(p_lims['min_s'], p_lims['max_s'])
        if normalize is True:
            plt.title('Normalized Bias in Fluxes for %s filter' % in_band)
        else:
            plt.title('Bias in Fluxes for %s filter' % in_band)
        plt.xlabel('5-sigma Depth (mags)')
        plt.ylabel('Observed Seeing (arcsec)')

        return fig

    def plot_bias_scatter(self, in_band='r', use_existing_fig=None,
                          normalize=False):
        """
        Plot the mean residuals of each visit individually as a function
        of depth and seeing.

        Parameters
        ----------
        in_band : str, default = 'r'
            The bandpass for the comparison. ('u', 'g', 'r', 'i', 'z', or 'y')
        use_existing_fig : matplotlib figure or None, default=None
            Can use an existing matplotlib figure or if set to None will
            return a new figure.
        normalize : Boolean, default=False
            If true, the plot values will be normalized by the standard
            deviation of the flux results in each bin.

        Returns
        -------
        fig : matplotlib figure
            Bias scatter plot with depth and seeing as axes.
        """
        p_lims, num_bins = self._set_stats_plot_limits(num_bins=0)
        in_band = str('lsst'+in_band)

        idx = np.where(self.flux_stats['bandpass'] == in_band)[0]

        if use_existing_fig is not None:
            fig = use_existing_fig
        else:
            fig = plt.figure()

        if normalize is True:
            c_vals = (self.flux_stats['mean_resid'][idx] /
                      np.sqrt((self.flux_stats['mean_sq_resid'][idx] -
                               self.flux_stats['mean_resid'][idx]**2.)))
            vmin = -1.
            vmax = 1.
        else:
            c_vals = (self.flux_stats['mean_resid'][idx])
            vmin = -0.3
            vmax = 0.3
        plt.scatter(self.flux_stats['depth'][idx],
                    self.flux_stats['seeing'][idx],
                    c=c_vals,
                    cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax, lw=0)
        plt.colorbar()
        plt.xlim(p_lims['min_d']-.05, p_lims['max_d']+.05)
        plt.ylim(p_lims['min_s']-.01, p_lims['max_s']+.01)
        if normalize is True:
            plt.title('Normalized Bias in Fluxes for %s bandpass' % in_band)
        else:
            plt.title('Bias in Fluxes for %s bandpass' % in_band)
        plt.xlabel('5-sigma Depth (mags)')
        plt.ylabel('Observed Seeing (arcsec)')

        return fig

    def plot_variance_map(self, in_band='r', with_bins=20.,
                          use_existing_fig=None):
        """
        Plot the mean of the squared residuals of each visit in a 2-d grid
        as a function of depth and seeing.

        Parameters
        ----------
        in_band : str, default = 'r'
            The bandpass for the comparison. ('u', 'g', 'r', 'i', 'z', or 'y')
        with_bins : float or length 2 list, default = 20.
            The number of bins to use to plot the depth and seeing values
            respectively.
        use_existing_fig : matplotlib figure or None, default=None
            Can use an existing matplotlib figure or if set to None will
            return a new figure.

        Returns
        -------
        fig : matplotlib figure
            Variance map with depth and seeing as axes.
        """
        p_lims, num_bins = self._set_stats_plot_limits(num_bins=with_bins)
        in_band = str('lsst'+in_band)

        d_grid, s_grid = np.meshgrid(np.linspace(p_lims['min_d']-.01,
                                                 p_lims['max_d']+.01,
                                                 num=num_bins[0]+1),
                                     np.linspace(p_lims['min_s']-.01,
                                                 p_lims['max_s']+.01,
                                                 num=num_bins[1]+1))
        p_lims['delta_d'] = d_grid[0, 1] - d_grid[0, 0]
        p_lims['delta_s'] = s_grid[1, 0] - s_grid[0, 0]
        p_val = np.zeros(np.shape(d_grid))
        for i in range(np.shape(d_grid)[0]):
            for j in range(np.shape(s_grid)[0]):
                depth_vals = self.flux_stats[((self.flux_stats['depth'] >=
                                               d_grid[i, j]) &
                                              (self.flux_stats['depth'] <
                                               (d_grid[i, j] +
                                                p_lims['delta_d'])) &
                                              (self.flux_stats['bandpass'] ==
                                               in_band))]
                grid_vals = depth_vals[((depth_vals['seeing'] >=
                                         s_grid[i, j]) &
                                        (depth_vals['seeing'] <
                                         (s_grid[i, j] + p_lims['delta_s'])))]

                if len(grid_vals) > 0:
                    p_val[i, j] = (np.average(grid_vals['mean_sq_resid'],
                                              weights=grid_vals[
                                                        'num_objects']) -
                                   np.average(grid_vals['mean_resid'],
                                              weights=grid_vals[
                                                        'num_objects'])**2.)

        if use_existing_fig is not None:
            fig = use_existing_fig
        else:
            fig = plt.figure()

        plt.pcolor(d_grid, s_grid, p_val, cmap=plt.cm.plasma,
                   vmin=0., vmax=0.45)
        plt.colorbar()
        plt.xlim(p_lims['min_d'], p_lims['max_d'])
        plt.ylim(p_lims['min_s'], p_lims['max_s'])
        plt.title('Flux Variance for %s bandpass' % in_band)
        plt.xlabel('5-sigma Depth (mags)')
        plt.ylabel('Observed Seeing (arcsec)')

        return fig

    def plot_variance_scatter(self, in_band='r', use_existing_fig=None):

        """
        Plot the mean of the squared residuals of each visit individually
        as a function of depth and seeing.

        Parameters
        ----------
        in_band : str, default = 'r'
            The bandpass for the comparison. ('u', 'g', 'r', 'i', 'z', or 'y')
        use_existing_fig : matplotlib figure or None, default=None
            Can use an existing matplotlib figure or if set to None will
            return a new figure.

        Returns
        -------
        fig : matplotlib figure
            Variance scatter plot with depth and seeing as axes.
        """

        p_lims, num_bins = self._set_stats_plot_limits(num_bins=0)
        in_band = str('lsst'+in_band)

        idx = np.where(self.flux_stats['bandpass'] == in_band)[0]

        if use_existing_fig is not None:
            fig = use_existing_fig
        else:
            fig = plt.figure()

        plt.scatter(self.flux_stats['depth'][idx],
                    self.flux_stats['seeing'][idx],
                    c=(self.flux_stats['mean_sq_resid'][idx] -
                       (self.flux_stats['mean_resid'][idx])**2),
                    cmap=plt.cm.plasma, vmin=0, vmax=0.45)

        plt.colorbar()
        plt.xlim(p_lims['min_d']-.05, p_lims['max_d']+.05)
        plt.ylim(p_lims['min_s']-.01, p_lims['max_s']+.01)
        plt.title('Flux Variance for %s bandpass' % in_band)
        plt.xlabel('5-sigma Depth (mags)')
        plt.ylabel('Observed Seeing (arcsec)')

        return fig

# =============================================================================


class BaseCurve(object):
    """
    A Base Class used to initialize the curve methods: LightCurve and
    SeeingCurve.

    ...

    Parameters
    ----------
    dbConn : dbInterface instance
        This is a connection to a science database with the LSST DM processed
        data products, e.g. a NERSC hosted pserv database.
        For an example see the 'depth_curve_example' notebook
        in the examples folder.
    fp_table_dir : str or None, default=None
        The location of a forced photometry fits table
    mjd_file : str or None, default=None
        The location of a CSV file with the visit numbers in the first column
        and the mjd values in the second.
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
            self.visit_map = {x: [] for x in self.bandpasses}
            self.lightcurve = None
            self.visit_mjd = {}

            # Associate mjd with visits (just a hack for now)
            # The whole ability to read from the data repo should probably
            # be updated to use the butler if we want to keep this ability
            mjd_list = np.genfromtxt(mjd_file, delimiter=',')
            for visit_date in mjd_list:
                self.visit_mjd[str(int(visit_date[0]))] = visit_date[1]

            for visit_dir in os.listdir(str(fp_table_dir+'/0/')):
                visit_band = visit_dir[-1]
                visit_num = visit_dir[1:-3]
                self.visit_map[visit_band].append(visit_num)


class LightCurve(BaseCurve):

    """
    Creates a pandas dataframe with the lightcurve information of an object.
    """

    def _build_lightcurve(self, lc_dict):
        """
        A wrapper around pd.read_dict to build a lightcurve stored in a
        pandas dataframe from a dict.
        """

        self.lightcurve = pd.DataFrame.from_dict(lc_dict)

    def build_lightcurve_from_fp_table(self, objid):
        """
        Assemble a light curve dataframe from available forced photometry
        in a data repo.

        Parameters
        ----------
        objid : int
            The object id for a source in the forced photometry tables.

        Attributes
        ----------
        lightcurve : pandas dataframe
            The lightcurve information stored in a pandas dataframe.
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
                    lightcurve['flux_error'].append(
                                         obj_data['base_PsfFlux_fluxSigma'][0])
                    lightcurve['zp'].append(25.0)
                    lightcurve['zpsys'].append('ab')
        self._build_lightcurve(lightcurve)

    def build_lightcurve_from_db(self, objid=None, ra_dec=None,
                                 tol=0.005):
        """
        Assemble a light curve from a pserv database by id or ra,dec location.

        Parameters
        ----------
        objid : int or None, default=None
            The object id for a source in the pserv database. If set to None
            then ra_dec must be specified.
        ra_dec : length 2 list or None, default=None
            The [ra, dec] location in degrees for objects in the pserv
            database. If set to None then objid must be specified. If there is
            more than one object that fits within the search radius set by
            tol then the possible object ids w/ ra, dec info
            will be printed out and the user can use objid to get
            lightcurves for the objects individually.
        tol : float, default=0.005
            The radius in degrees to search in ra, dec locations for
            objects in the pserv database.

        Attributes
        ----------
        lightcurve : pandas dataframe
            The lightcurve information stored in a pandas dataframe.
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
            obj_info = self.dbConn_lc.objectFromRaDec(ra_dec[0], ra_dec[1],
                                                      tol)
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
        lightcurve['mag_error'] = 2.5*np.log10(1 + (fs_info['psf_flux_err'] /
                                                    fs_info['psf_flux']))
        lightcurve['zp'] = [25.0]*num_results  # TEMP
        lightcurve['zpsys'] = ['ab']*num_results  # TEMP

        self._build_lightcurve(lightcurve)

    def visualize_lightcurve(self, using='flux', include_errors=True,
                             use_existing_fig=None):
        """
        Make a simple light curve plot.

        Adapted from sncosmo.plot_lc source code.

        Parameters
        ----------
        using : str, ('flux' or 'mag'), default='flux'
            Plot the flux or the magnitude of the object over time.
        include_errors : boolean, default=True
            Set to true to include error bars on the plot.
        use_existing_fig : matplotlib figure or None, default=None
            Specify an existing matplot to use or set to None to make one.

        Returns
        -------
        fig : matplotlib figure
            Plot with object flux or magnitude over time.
        """
        if self.lightcurve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        n_subplot = len(self.filter_list)
        n_col = 2
        n_row = (n_subplot - 1) // n_col + 1
        if use_existing_fig is None:
            fig = plt.figure(figsize=(4. * n_col, 3. * n_row))
        else:
            fig = use_existing_fig

        color = ['b', 'g', 'y', 'orange', 'r', 'k']

        plot_num = 1
        for filt in self.filter_list:
            fig.add_subplot(n_row, n_col, plot_num)
            filt_name = str('lsst' + filt)
            plt.title(filt_name)
            dc_bp = self.lightcurve['bandpass']
            filt_mjd = self.lightcurve['mjd'][dc_bp == filt_name].values
            using_mjd = self.lightcurve[using][dc_bp == filt_name].values
            if include_errors is True:
                err_str = str(using+'_error')
                using_error_mjd = self.lightcurve[err_str][dc_bp ==
                                                           filt_name].values
                plt.errorbar(filt_mjd, using_mjd, yerr=using_error_mjd,
                             ls='None', marker='.', ms=3, c=color[plot_num-1],
                             capsize=4)
            else:
                plt.scatter(filt_mjd, using_mjd, c='purple', marker='+')
            plt.locator_params(axis='x', nbins=5)
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

# =============================================================================


class SeeingCurve(BaseCurve):
    """
    Store the seeing information from all visits in the survey in a pandas
    dataframe.
    """

    def _build_seeing_curve(self, sc_dict):
        """
        A wrapper around pd.read_dict to build a seeing curve stored in a
        pandas dataframe from a dict.
        """

        self.seeing_curve = pd.DataFrame.from_dict(sc_dict)

    def visualize_seeing_curve(self):

        if self.seeing_curve is None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        n_col = 2
        n_row = len(self.filter_list)
        fig = plt.figure(figsize=(4. * n_col, 3. * n_row))

        color = ['b', 'g', 'y', 'orange', 'r', 'k']

        plot_num = 1
        for filt in self.filter_list:

            fig.add_subplot(n_row, n_col, plot_num)
            filt_name = str('lsst' + filt)
            plt.title(filt_name)

            sc_bp = self.seeing_curve['bandpass']
            filt_seeing = self.seeing_curve['seeing'][sc_bp == filt_name]
            filt_mjd = self.seeing_curve['mjd'][sc_bp == filt_name]
            plt.scatter(filt_mjd, filt_seeing,
                        c=color[(plot_num - 1)//2], marker='+')
            plt.locator_params(axis='x', nbins=5)
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
