# ==============================================================================
# License info here?
# ==============================================================================

import monitor
import os
import lsst.daf.persistence as dp
import pandas as pd
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

    def read(self,photfile):
        print("No IO enabled yet.")
        return

    def run(self,algorithm=None):
        if algorithm is None:
            print("No algorithms coded yet.")
        return

class LightCurve(object):

    '''
    Fetch ForcedSource data with Butler and package for use with accompanying visualization
    routines.
    '''

    def __init__(self, fp_table_dir=None, bandpasses=None, visitLists=None, mjdFile=None):

        self.fp_table_dir = fp_table_dir
        self.visitMap = {}
        self.bandpasses = bandpasses
        self.lightcurve = None
        self.visitMjd={}

        #Associate mjd with visits (just a hack for now, need to think about this more)
        mjdList = np.genfromtxt(mjdFile,delimiter=',')
        for visitDate in mjdList:
            self.visitMjd[str(int(visitDate[0]))] = visitDate[1]

        for bandpass, visitList in zip(self.bandpasses, visitLists):
            #Associate correct visits with bandpass
            self.visitMap[bandpass] = visitList
            #Register required lsst bandpass in sncosmo registry
            bandpassFile = os.path.join(str(getPackageDir('throughputs') + '/baseline/total_' +
                                            bandpass + '.dat'))
            bandpassInfo = np.genfromtxt(bandpassFile, names=['wavelen', 'transmission'])
            band = sncosmo.Bandpass(bandpassInfo['wavelen'], bandpassInfo['transmission'],
                                    name=str('lsst' + bandpass))
            sncosmo.registry.register(band)


    def build_lightcurve(self, objid):

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
            for visit in self.visitMap[bandpass]:
                hdulist = fits.open(str(self.fp_table_dir + '/0/v' + str(visit) + '-fr/R22/S11.fits'))
                objData = hdulist[1].data[np.where(hdulist[1].data['objectId']==objid)]
                if len(objData) > 0:
                    lightcurve['bandpass'].append(str('lsst' + bandpass))
                    lightcurve['mjd'].append(self.visitMjd[str(visit)])
                    lightcurve['ra'].append(objData['coord_ra'][0])
                    lightcurve['dec'].append(objData['coord_dec'][0])
                    lightcurve['flux'].append(objData['base_PsfFlux_flux'][0])
                    lightcurve['flux_error'].append(objData['base_PsfFlux_fluxSigma'][0])
                    lightcurve['zp'].append(25.0)
                    lightcurve['zpsys'].append('ab')
        self.lightcurve = lightcurve

    def visualize_lightcurve(self):

        if self.lightcurve == None:
            raise ValueError('No lightcurve yet. Use build_lightcurve first.')

        lcTable = Table(data=self.lightcurve)
        fig = sncosmo.plot_lc(lcTable)

        return fig





# ==============================================================================
# Need better tests...

if __name__ == '__main__':

    newMonitor = monitor.Monitor()

# ==============================================================================
