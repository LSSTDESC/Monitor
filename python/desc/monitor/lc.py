"""
Class for representing light curve
"""
import pandas as pd

class LC(object):

    def __init__(self, mjds, zps, zpsys, flux, flxuerrors, bands, sourceIDs):

        self.mjds = mjds
        self.zps = zps
        self.zpsys = zpsys
        self.flux = flux
        self.fluxerr = fluxerr
        self.data = pd.DataFrame({'time': self.mjds,
                                  'flux': self.flux,
                                  'fluxerr': self.fluxerr,
                                  'zp': self.zps
                                  'zpsys': self.zpsys})


    @classmethod
    def fromSourceID(cls, )


