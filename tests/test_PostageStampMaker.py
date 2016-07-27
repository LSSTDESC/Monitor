"""
Test code for PostageStampMaker module.
"""
from __future__ import absolute_import, division
import os
import unittest
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.utils
from desc.monitor import PostageStampMaker, create_postage_stamps,\
    convert_image_to_hdu

class PostageStampTestCase(unittest.TestCase):
    "TestCase class for PostageStampMaker module."
    def setUp(self):
        """
        Set up each test with a new PostageStampMaker object using the
        test FITS file, coordinates of the image center, and a postage
        stamp size of 10 arcsec.
        """
        self.expfile = os.path.join(lsst.utils.getPackageDir('monitor'),
                                    'tests', 'small_CoaddTempExp.fits.gz')
        self.stamp_maker = PostageStampMaker(self.expfile)
        center_coord = self.stamp_maker.center_coord(self.stamp_maker.exposure)
        self.ra = center_coord.getLongitude().asDegrees()
        self.dec = center_coord.getLatitude().asDegrees()
        self.size = 10

    def tearDown(self):
        "Clean up"
        del self.stamp_maker

    def test_bbox_generation(self):
        "Test the make bounding box function."
        bbox = self.stamp_maker.make_bbox(self.ra, self.dec, self.size)
        self.assertEqual((bbox.getWidth(), bbox.getHeight()), (30, 30))

    def test_stamp_size(self):
        "Test that the postage stamp has the expected size."
        my_exposure = self.stamp_maker.create(self.ra, self.dec, self.size)
        my_imarr = my_exposure.getMaskedImage().getImage().getArray()
        self.assertEqual(my_imarr.shape, (30, 30))

    def test_stamp_central_pixel_value(self):
        """
        Test that the center pixel of both the postage stamp and original
        image have the same value.
        """
        my_exposure = self.stamp_maker.create(self.ra, self.dec, self.size)
        my_imarr = my_exposure.getMaskedImage().getImage().getArray()
        ref_exposure = self.stamp_maker.exposure
        ref_imarr = ref_exposure.getMaskedImage().getImage().getArray()
        self.assertEqual(my_imarr[int(my_imarr.shape[0]/2)][int(my_imarr.shape[1]/2)],
                         ref_imarr[int(ref_imarr.shape[0]/2)][int(ref_imarr.shape[1]/2)])

    def test_stamp_centers_match(self):
        """
        Test that coordinates of the centers of the postage stamp and
        the original image are the same.
        """
        my_exposure = self.stamp_maker.create(self.ra, self.dec, self.size)
        self.assertEqual(self.stamp_maker.center_coord(my_exposure),
                         self.stamp_maker.center_coord(self.stamp_maker.exposure))

    def test_create_array_of_stamps(self):
        """
        Test that creates a sequence of stamps from the Exposure object
        given a sequence of (ra, dec, size) tuples and tests that the
        centers of the stamps are at the expected locations.
        """
        pix_scale = self.stamp_maker.exposure.getWcs().pixelScale().asDegrees()
        npix = 10
        stamp_specs = [(self.ra, self.dec - npix*pix_scale, self.size),
                       (self.ra, self.dec, self.size),
                       (self.ra, self.dec + npix*pix_scale, self.size)]
        my_stamps = self.stamp_maker.create_stamps(stamp_specs)
        self.assertEqual(len(my_stamps), len(stamp_specs))
        for i, stamp in enumerate(my_stamps):
            coord = self.stamp_maker.center_coord(stamp)
            self.assertAlmostEqual(coord.getLongitude().asDegrees(),
                                   stamp_specs[i][0])
            self.assertAlmostEqual(coord.getLatitude().asDegrees(),
                                   stamp_specs[i][1])

    def test_create_sequence_function(self):
        """
        For a given (ra, dec, size), test the create_postage_stamps
        function which returns a sequence of stamps for that sky
        region, extracted from a sequence of input Exposure FITS
        files.
        """
        fits_files = [self.expfile]*3
        my_stamps = create_postage_stamps(self.ra, self.dec, self.size,
                                          fits_files)
        self.assertEqual(len(my_stamps), len(fits_files))
        for stamp in my_stamps:
            self.assertTrue(isinstance(stamp, afwImage.ExposureF))

    def test_hdu_from_exposure(self):
        """
        Test the conversion of the image from Exposure object directly
        to an HDU.
        """
        my_stamp = self.stamp_maker.create(self.ra, self.dec, self.size)
        my_hdu = convert_image_to_hdu(my_stamp)

        # Write the Exposure to a FITS file to test the image data and
        # WCS keywords in the in-memory HDU against the versions
        # written by the Stack.
        stamp_file = 'temp_stamp.fits'
        my_stamp.writeFits(stamp_file)
        stack_hdu = fits.open(stamp_file)[1]

        for value1, value2 in zip(my_hdu.data.flatten(),
                                  stack_hdu.data.flatten()):
            self.assertEqual(value1, value2)
        for keyword in my_hdu.header.keys():
            self.assertEqual(my_hdu.header[keyword], stack_hdu.header[keyword])
        os.remove(stamp_file)

if __name__ == '__main__':
    unittest.main()
