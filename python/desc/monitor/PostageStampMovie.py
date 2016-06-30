"""
Use matplotlib.animation tools to make a movie of Twinkles data.
"""
from __future__ import print_function, absolute_import, division
import os
import sys
import copy
import pickle
from collections import OrderedDict
import numpy as np
import astropy.time
import matplotlib.pyplot as plt
from matplotlib import animation
from desc.pserv import DbConnection
from .PostageStampMaker import PostageStampMaker, convert_image_to_hdu
from .Display import image_norm, render_fits_image
plt.ion()

__all__ = ['PostageStampMovie', 'Level2DataService']

class Level2DataService(object):
    """Provide access to the Twinkles Level 2 data from a MySQL/pserv
    database.

    Parameters
    ----------
    repo : str
        The output repository from the Level 2 analysis.
    db_info : dict or None, optional
        Connection information for the Level 2 database.  It should
        have the form

        {'database': <database name>,
         'host': <database server host address>}

        If None, then the default info for using the Twinkles Level 2
        database at NERSC is used:

        {'database': 'DESC_Twinkles_Level_2',
         'host': 'scidb1.nersc.gov'}

    Attributes
    ----------
    repo : str
        The output repository from the Level 2 analysis.
    conn : desc.pserv.DbConnection
        The connection object that executes the queries to the Level 2 db.

    Notes
    -----
    Database credentials and connection info (e.g., the port to use)
    are accessed via lsst.daf.persistence.DbAuth from the user's
    ~/.lsst/db-auth.paf file.

    """
    def __init__(self, repo, db_info=None):
        """
        Class constructor.
        """
        self.repo = repo
        if db_info is None:
            db_info = dict(database='DESC_Twinkles_Level_2',
                           host='scidb1.nersc.gov')
        self.conn = DbConnection(**db_info)

    def get_pixel_data(self, objectId, band, size=10, pickle_file=None):
        """Get the pixel data for a cutout region for all visits.

        Use PostageStampMaker to extract the pixel data cutout for the
        specified objectId and band for each visit in the repository.

        Parameters
        ----------
        objectId : int
           The id of the object in the Level2 Object db table.
        band : str
           The LSST band (u, g, r, i, z, or y).
        size : int, optional
           The size of the cutout region in arcsec.
           Cutouts of size x size centered on the RA, Dec of the
           object will be extracted.
        pickle_file : str, optional
           The name of the pickle file to contain the extracted pixel
           data.  If it is None, then a default name of the form
           'pixel_data_%(objectId)i_%(band)s.pkl' will be used.  If
           that file exists, the pixel data will be loaded from it.
           If it does not exist, the pixel data will be extracted from
           the deepCoadd tempExp frames.

        Returns
        -------
        pixel_data : OrderedDict
            A dictionary of 2D numpy arrays containing the pixel values.
            The keys are the visitId numbers from the CcdVisit table.
        pickle_file : str
            The name of the pickle file.
        """
        if pickle_file is None:
            pickle_file = 'pixel_data_%(objectId)i_%(band)s.pkl' % locals()
        if os.path.isfile(pickle_file):
            with open(pickle_file, 'r') as input_:
                pixel_data = pickle.load(input_)
            return pixel_data, pickle_file
        pixel_data = OrderedDict()
        visits = self.get_visits(band)
        ra, dec = self.get_coords(objectId)
        for visit in visits:
            print("working on visit", visit)
            sys.stdout.flush()
            exp_file = self._get_warped_visit_frame(band, visit)
            maker = PostageStampMaker(exp_file)
            stamp = maker.create(ra, dec, size)
            pixel_data[visit] = \
                copy.deepcopy(stamp.getMaskedImage().getImage().getArray())
        with open(pickle_file, 'w') as output:
            pickle.dump(pixel_data, output)
        return pixel_data, pickle_file

    def _get_warped_visit_frame(self, band, visit):
        """Get the deepCoadd tempExp frame.

        Parameters
        ----------
        band : str
            The LSST band (u, g, r, i, z, or y).
        visit : int
            The visitId of the visit from the CcdVisit table.

        Returns
        -------
        str
            File path to the FITS file containing the warped exposure for
            the specified band and visit.

        Notes
        -----
        This should be re-implemented to use the Data Butler.
        """
        return os.path.join(self.repo, 'deepCoadd', band, '0/0,0tempExp',
                            'v%(visit)i-f%(band)s.fits' % locals())

    def get_light_curve(self, objectId, band):
        """ Get the light curve for the requested objectId and band.

        Parameters
        ----------
        objectId : int
            The id of the object as specified in the Level 2 Object table.
        band : str
            The LSST band (u, g, r, i, z, or y).

        Returns
        -------
        dict
            A dictionary of numpy arrays containing the light curve data.
            The format is
            {'mjd': <Visit obs times in MJD>,
             'flux': <Fluxes in nmgy>,
             'fluxerr': <Flux errors in nmgy>}

        Notes
        -----
        This method queries the CcdVisit and ForcedSource tables for
        the CcdVisit.obsStart, ForcedSource.psFlux, and
        ForcedSource.psFlux_Sigma values.
        """
        query = """select cv.obsStart, fs.psFlux, fs.psFlux_Sigma from
               CcdVisit cv join ForcedSource fs on cv.ccdVisitId=fs.ccdVisitId
               join Object obj on fs.objectId=obj.objectId where
               cv.filterName='%(band)s' and fs.objectId=%(objectId)i
               order by cv.obsStart asc""" % locals()
        rows = self.conn.apply(query, lambda curs: np.array([x for x in curs]))
        obsStart, flux, fluxerr = (np.array(col) for col in rows.transpose())
        mjd = astropy.time.Time(obsStart).mjd
        return dict(mjd=mjd, flux=flux, fluxerr=fluxerr)

    def get_visits(self, band):
        """Get the visitIds corresponding to the specified band.

        Parameters
        ----------
        band : str
            The LSST band (u, g, r, i, z, or y).

        Returns
        -------
        list
            A list of visitIds from the CcdVisit table.
        """
        return self.conn.apply('''select visitId from CcdVisit where
                                  filterName='%(band)s' order by visitId'''
                               % locals(),
                               lambda curs: [x[0] for x in curs])

    def get_coords(self, objectId):
        """Get the RA, Dec of the requested object from the Object db table.

        Parameters
        ----------
        objectId : int
            The id of the object as specified in the Level 2 Object table.

        Returns
        -------
        tuple
            The objects RA, Dec in degrees.
        """
        return self.conn.apply('''select psRA, psDecl from Object where
                                  objectId=%(objectId)i''' % locals(),
                               lambda curs: tuple([x for x in curs][0]))

class PostageStampMovie(object):
    """Produce an animation of a cutout around a coadd object.

    The forced source light curve is also displayed and cursor events
    on the light curve plot can be used to control the animation.

    Parameters
    ----------
    objectId : int
        The id of the object in the Level2 Object db table.
    band : str
        The LSST band (u, g, r, i, z, or y).
    size : int, optional
        The size of the cutout region in arcsec.
        Cutouts of size x size centered on the RA, Dec of the
        object will be extracted.
    pickle_file : str, optional
        The name of the pickle file to contain the extracted pixel
        data.  If it is None, then a default name of the form
        'pixel_data_%(objectId)i_%(band)s.pkl' will be used.  If
        that file exists, the pixel data will be loaded from it.
        If it does not exist, the pixel data will be extracted from
        the deepCoadd tempExp frames.
    figsize : (float, float), optional
        The size of the figure produced by matplotlib.
    scaling_factor : float
        The scaling factor to apply to the image normalization applied to
        each visit frame.  The image normalization is inferred from the
        coadd frame.  See the documentation for desc.monitor.image_norm.

    Attributes
    ----------
    objectId : int
        The id of the object in the Level2 Object db table.
    band : str
        The LSST band (u, g, r, i, z, or y).
    pixel_data : OrderedDict
        A dictionary of 2D numpy arrays containing the pixel values.
        The keys are the visitId numbers from the CcdVisit table.
    pickle_file : str
        The name of the pickle file to contain the extracted pixel data.
    fig : matplotlib.figure.Figure
        The Figure object containing the cutout and light curve plot.
    image : matplotlib.image.AxesImage
        The AxesImage object displaying the cutout animation.
    light_curve : dict
        A dictionary of numpy arrays containing the light curve data.
        The format is
        {'mjd': <Visit obs times in MJD>,
         'flux': <Fluxes in nmgy>,
         'fluxerr': <Flux errors in nmgy>}
    """
    def __init__(self, objectId, band, l2_service, size=10,
                 pickle_file=None, figsize=(6, 10), scaling_factor=50):
        """
        Constructor to create the figure with the animated cutout
        and lightcurve.
        """
        self.objectId = objectId
        self.band = band
        self.pixel_data, self.pickle_file = \
            l2_service.get_pixel_data(objectId, band, size=size,
                                      pickle_file=pickle_file)
        self._display_figure(l2_service, size, figsize, scaling_factor)
        self._set_animation_attributes()
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def _display_figure(self, l2_service, size, figsize, scaling_factor):
        """Display the image and light curve.

        Parameters
        ----------
        l2_service : Level2DataService object
            The interface to the Level 2 data.
        size : int
            The linear size of the square cutout in arcsec.
        figsize : (float, float)
            The size of the matplotlib.figure.Figure containing the animation.
        scaling_factor : float
            The scaling factor to apply to the image normalization
            applied to each visit frame.  The image normalization is
            inferred from the coadd frame.  See the documentation for
            desc.monitor.image_norm.
        """
        objectId = self.objectId
        band = self.band

        # Create the coadd postage stamp and use as the initial image.
        coadd = PostageStampMaker(os.path.join(l2_service.repo, 'deepCoadd',
                                               band, '0/0,0.fits'))
        ra, dec = l2_service.get_coords(objectId)

        # Use the coadd to set the image normalization
        stamp = coadd.create(ra, dec, size)
        hdu = convert_image_to_hdu(stamp)
        norm = image_norm(hdu.data*scaling_factor)

        plt.rcParams['figure.figsize'] = figsize
        self.fig = plt.figure()
        self.fig.suptitle('objectId: %(objectId)i\n%(band)s band' % locals())
        axes, self.image = render_fits_image(hdu, norm=norm, fig=self.fig,
                                             subplot=211)[1:3]
        axis_range = plt.axis()
        axes.scatter([ra], [dec], transform=axes.get_transform('icrs'),
                     color='red', alpha=0.8)
        plt.axis(axis_range)   # restore axis range to fit image
        self.fig.add_subplot(212)
        self.light_curve = l2_service.get_light_curve(objectId, band)
        plt.errorbar(self.light_curve['mjd'], self.light_curve['flux'],
                     yerr=self.light_curve['fluxerr'], fmt='.')
        self._yrange = plt.axis()[2:]
        self._current_point = plt.plot([self.light_curve['mjd'][0]],
                                       [self.light_curve['flux'][0]],
                                       marker='o', color='red')
        self._current_point.extend(plt.plot([self.light_curve['mjd'][0],
                                             self.light_curve['mjd'][0]],
                                            self._yrange, 'k:'))
        plt.xlabel('MJD')
        plt.ylabel('flux (nmgy)')

    def _set_animation_attributes(self):
        "Set the initial values of the attributes to control the animation."
        self._pause = False
        self._update = False
        self._num = 0
        self._nmax = len(self.light_curve['mjd'])

    def run(self, interval=200):
        """Run the animation with a time between frames of interval.

        The returned FuncAnimation object must exist in the top-level
        context.

        Parameters
        ----------
        interval : int
            The interval time between frames in msec

        Return
        ------
        matplotlib.animation.FuncAnimation
        """
        return animation.FuncAnimation(self.fig, self, frames=self.index,
                                       interval=interval)

    def index(self):
        """Generator that returns the frame number to display.

        Yields
        ------
        int
            The frame number to display.
        """
        while True:
            yield self._num
            if not self._pause or self._update:
                self._num += 1
            if self._num >= self._nmax or self._num < 0:
                self._num = 0

    def __call__(self, i):
        """Call-back to set the data in the image for the ith frame.

        Parameters
        ----------
        i : int
            The frame number to display.

        Returns
        -------
        [matplotlib.image.AxesImage,
         [matplotlib.axes.Axes, matplotlib.axes.Axes]]
        """
        mjd = self.light_curve['mjd']
        flux = self.light_curve['flux']
        if not self._pause or self._update:
            self.image.set_data(self.pixel_data.values()[i])
            self._current_point[0].set_data([mjd[i]], [flux[i]])
            self._current_point[1].set_data([mjd[i], mjd[i]], self._yrange)
            self._update = False
        return [self.image, self._current_point]

    def on_click(self, event):
        """Method to set mouse event info.

        Parameters
        ----------
        event : matplotlib mouse event
        """
        if event.button == 3:
            try:
                dt = np.abs(self.light_curve['mjd'] - event.xdata)
                self._num = np.where(dt == min(dt))[0][0] - 1
                self._update = True
            except IndexError:
                pass
        if event.button == 2:
            self._pause ^= True

if __name__ == '__main__':
    repo = '/nfs/farm/g/desc/u1/users/jchiang/desc_projects/twinkles/Run1.1/output'
    band = 'r'
    interval = 200
    db_info = dict(host='ki-sr01.slac.stanford.edu',
                   port='3307',
                   database='jc_desc')
#    db_info = dict(host='127.0.0.1',
#                   port='53306',
#                   database='DESC_Twinkles_Level_2')
    level2_service = Level2DataService(repo, db_info=db_info)
    movies = []
    for objectId in (6931, 50302, 52429)[:1]:
        ps = PostageStampMovie(objectId, band, level2_service)
        movies.append(ps.run(interval))

#    objectId = 48237
#    for band in 'ur':
#        ps = PostageStampMovie(objectId, band, level2_service)
#        movies.append(ps.run(interval))
