from __future__ import absolute_import, print_function
import numpy as np
import sys
from PyQt4 import QtGui, QtCore

# Load in the light curve of a given supernova

mjds, mags = np.loadtxt("./SN_opsim.csv", dtype='string', delimiter=',',
                        usecols=(1, 4), unpack=True)
mjds = mjds.astype("double")
mags = mags.astype("double")
mjds = mjds[~np.isnan(mags)]
mjds = mjds.astype("int")


class Example(QtGui.QWidget):

    def __init__(self):
        super(Example, self).__init__()

        self.initUI()

    def initUI(self):

        # ------------------------------------------
        # Initialize the sliderbar, and connect the value of the slider and the
        # images.
        #
        sld = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        sld.setFocusPolicy(QtCore.Qt.NoFocus)
        sld.setGeometry(40, 700, 340, 40)
        sld.setRange(0, 239)
        sld.valueChanged[int].connect(self.changeValue)

        # ------------------------------------------
        # Initialize the lcd number for expMJDs, and connect the value of the
        # slider and the expMJD value.
        #
        self.lcd = QtGui.QLCDNumber(self)
        sld.valueChanged[int].connect(self.showMjds)

        # ------------------------------------------
        # Initialize the pictures of the light curve
        #
        self.label1 = QtGui.QLabel(self)
        pixmap0 = QtGui.QPixmap('./pngs/pngs0000000000_sn.png')
        pixmap = pixmap0.scaled(340, 340, QtCore.Qt.KeepAspectRatio)
        self.label1.setPixmap(pixmap)
        self.label1.setGeometry(40, 0, 340, 340)

        # ------------------------------------------
        # Initialize the pictures of the images of the supernova
        #
        self.label2 = QtGui.QLabel(self)
        pixmap0 = QtGui.QPixmap('./pngs/pngs0000000000_sn.png')
        pixmap = pixmap0.scaled(340, 340, QtCore.Qt.KeepAspectRatio)
        self.label2.setPixmap(pixmap)
        self.label2.setGeometry(40, 320, 340, 340)

        # ------------------------------------------
        # Initialize the presentation of expMJDs
        #
        self.lcd.display(mjds[0])
        self.lcd.setGeometry(200, 660, 100, 40)
        self.label0 = QtGui.QLabel(self)
        self.label0.setText('expMJDs:')
        self.label0.setFont(QtGui.QFont('SansSerif', 28))
        self.label0.setGeometry(50, 660, 200, 40)

        # ------------------------------------------
        # Initialize the window of this visualization
        #
        self.setGeometry(500, 500, 420, 768)
        self.setWindowTitle('Light Curve of a Supernova')
        self.show()

    def changeValue(self, value):
        # ------------------------------------------
        # Pictures of the light curve
        #
        filename = './pngs/pngs'+'{:0>10}'.format(str(value))+'_sn.png'
        pixmap0 = QtGui.QPixmap(filename)
        pixmap = pixmap0.scaled(340, 340,QtCore.Qt.KeepAspectRatio)
        self.label1.setPixmap(pixmap)

        # ------------------------------------------
        # Pictures of the images of supernova
        #
        filename = './pngs/pngs'+'{:0>10}'.format(str(value))+'_sn.png'
        pixmap0 = QtGui.QPixmap(filename)
        pixmap = pixmap0.scaled(340, 340, QtCore.Qt.KeepAspectRatio)
        self.label2.setPixmap(pixmap)

    def showMjds(self, value):
        self.lcd.display(mjds[value])


def main():

    app = QtGui.QApplication(sys.argv)
    # ex = Example()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
