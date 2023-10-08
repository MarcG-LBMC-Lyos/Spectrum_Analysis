# SpectrumAnalysis is a program of automatic analysis of bone FTIRM spectrum.
# Copyright (C) 2017  GARDEGARONT Marc
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from tkinter.filedialog import asksaveasfilename


def writeDefaultConfig(filename, writeCurrent=False, **kwargs):
    """ Writes the default config file """
    if not writeCurrent:
        MMABOUNDS = [1706, 1826]

        v1v3PO4 = [910, 1184]
        guess_v1v3PO4 = [0.0001, 5, 962, 0.0001, 5, 1030, 0.0001, 5, 1060,
                         0.0001, 5, 1082, 0.0001, 5, 1110]
        bounds_v1v3PO4 = [[0.0001, 0.0001, 952, 0.0001, 0.0001, 1022, 0.0001, 0.0001,
                           1045, 0.0001, 0.0001, 1072, 0.0001, 0.0001, 1100],
                          [np.Inf, 35.000, 972, np.Inf, 35.000, 1038, np.Inf, 25.000,
                           1075, np.Inf, 25.000, 1092, np.Inf, 35.000, 1120]]

        v4PO4 = [500, 650]
        guess_v4PO4 = [0.0001, 5, 552, 0.0001, 5, 563, 0.0001, 5, 577,
                       0.0001, 5, 604]
        bounds_v4PO4 = [[0.0001, 0.0001, 546, 0.0001, 0.0001, 553, 0.0001, 0.0001,
                         567, 0.0001, 0.0001, 594],
                        [np.Inf, np.Inf, 558, np.Inf, 20.000, 573, np.Inf, 20.000,
                         587, np.Inf, np.Inf, 614]]

        CO3 = [862, 894]
        guess_CO3 = [0.0001, 5, 866, 0.0001, 5, 873, 0.0001, 5, 880]
        bounds_CO3 = [[0.0001, 0.0001, 864, 0.0001, 0.0001, 870, 0.0001, 0.0001, 876],
                      [np.Inf, np.Inf, 868, np.Inf, np.Inf, 876, np.Inf, np.Inf, 884]]

        amide = [1300, 1817]
        guess_amide = [0.0001, 5, 1340, 0.0001, 5, 1420, 0.0001, 5, 1450, 0.0001, 5, 1550, 0.0001, 5, 1600, 0.0001, 5, 1633, 0.0001, 5, 1660, 0.0001, 5, 1690]
        bounds_amide = [[0.0001, 0.0001, 1330, 0.0001, 0.0001, 1410, 0.0001, 0.0001, 1440, 0.0001, 0.0001, 1540, 0.0001, 0.0001, 1590, 0.0001, 0.0001, 1630, 0.0001, 0.0001, 1658, 0.0001, 0.0001, 1689],
                        [np.Inf, 35.000, 1350, np.Inf, 35.000, 1430, np.Inf, 35.000,
                         1460, np.Inf,
                         35.000, 1560, np.Inf, 35.000, 1610, np.Inf, 25.000, 1636,
                         np.Inf, 35.000, 1662, np.Inf, 35.000, 1691]]

    else:
        MMABOUNDS = kwargs.get('mmabounds', None)

        v1v3PO4 = kwargs.get('v1v3PO4', None)
        guess_v1v3PO4 = kwargs.get('guess_v1v3PO4', None)
        bounds_v1v3PO4 = kwargs.get('bounds_v1v3PO4', None)

        v4PO4 = kwargs.get('v4PO4', None)
        guess_v4PO4 = kwargs.get('guess_v4PO4', None)
        bounds_v4PO4 = kwargs.get('bounds_v4PO4', None)

        CO3 = kwargs.get('CO3', None)
        guess_CO3 = kwargs.get('guess_CO3', None)
        bounds_CO3 = kwargs.get('bounds_CO3', None)

        amide = kwargs.get('amide', None)
        guess_amide = kwargs.get('guess_amide', None)
        bounds_amide = kwargs.get('bounds_amide', None)

    with open(filename, 'w+') as conffile:
        # Comments
        conffile.write(
            '# Configuration file for automatic ' +
            'deconvolution. function is the function used to deconvolve the ' +
            'signal (either gaussian or lorentzian (voigt is usable but not right)). ' +
            'Xbound_sub are the boundaries used ' +
            'to detect and remove the X spectrum in the total ' +
            'spectrum (e.g. approximate boundaries of the highest MMA peak).' +
            'Possibilities to subtract more than one spectrum.'
            'The following structure is as this : consider X the' +
            'the name of a Region Of Interest in the spectrum. ' +
            'Xbound are the boundaries of this ROI. guess_height_X ' +
            'are the initial guess for the height of each peak in ' +
            'X. guess_halffwhm is the same but for the half FWHM. ' +
            'And guess_wavenb also, but this time for the position ' +
            'of the peak (wave number). boundsinf_*_X have the same ' +
            'meaning but for the inferior boundary (for the ' +
            'optimization process / deconvolution). And finally, ' +
            'boundssup_*_X also, but for the superior boundary.\n')
        conffile.write(
            '# WARNING : IF YOU WANT TO PUT A NULL VALUE, DONT PUT ' +
            '0 BUT PREFER A NEAR ZERO VALUE LIKE 0.000001. IF YOU ' +
            'WANT TO PUT AN INFINITY VALUE, JUST WRITE inf.\n')

        # Function used to deconvolve the spectrum
        conffile.write('function gaussian\n')

        # Limit of absorbance for the spectrometer captor
        conffile.write('limit 1.5\n')

        # MMA peak boundaries
        conffile.write('MMAbound_sub {0} {1}\n'.format(MMABOUNDS[0],
                                                       MMABOUNDS[1]))
        conffile.write('\n')

        # v1v3PO4 area
        conffile.write('v1v3PO4bound {0} {1}\n'.format(v1v3PO4[0], v1v3PO4[1]))
        conffile.write('guess_height_v1v3PO4')
        [conffile.write(' {}'.format(guess_v1v3PO4[i * 3]))
            for i in range(len(guess_v1v3PO4) // 3)]
        conffile.write('\n')
        conffile.write('guess_halffwhm_v1v3PO4')
        [conffile.write(' {}'.format(guess_v1v3PO4[i * 3 + 1]))
            for i in range(len(guess_v1v3PO4) // 3)]
        conffile.write('\n')
        conffile.write('guess_wavenb_v1v3PO4')
        [conffile.write(' {}'.format(guess_v1v3PO4[i * 3 + 2]))
            for i in range(len(guess_v1v3PO4) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_height_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[0][i * 3]))
            for i in range(len(bounds_v1v3PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_height_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[1][i * 3]))
            for i in range(len(bounds_v1v3PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_halffwhm_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[0][i * 3 + 1]))
            for i in range(len(bounds_v1v3PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_halffwhm_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[1][i * 3 + 1]))
            for i in range(len(bounds_v1v3PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_wavenb_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[0][i * 3 + 2]))
            for i in range(len(bounds_v1v3PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_wavenb_v1v3PO4')
        [conffile.write(' {}'.format(bounds_v1v3PO4[1][i * 3 + 2]))
            for i in range(len(bounds_v1v3PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('\n')

        # v4PO4 area
        conffile.write('v4PO4bound {0} {1}\n'.format(v4PO4[0], v4PO4[1]))
        conffile.write('guess_height_v4PO4')
        [conffile.write(' {}'.format(guess_v4PO4[i * 3]))
            for i in range(len(guess_v4PO4) // 3)]
        conffile.write('\n')
        conffile.write('guess_halffwhm_v4PO4')
        [conffile.write(' {}'.format(guess_v4PO4[i * 3 + 1]))
            for i in range(len(guess_v4PO4) // 3)]
        conffile.write('\n')
        conffile.write('guess_wavenb_v4PO4')
        [conffile.write(' {}'.format(guess_v4PO4[i * 3 + 2]))
            for i in range(len(guess_v4PO4) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_height_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[0][i * 3]))
            for i in range(len(bounds_v4PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_height_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[1][i * 3]))
            for i in range(len(bounds_v4PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_halffwhm_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[0][i * 3 + 1]))
            for i in range(len(bounds_v4PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_halffwhm_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[1][i * 3 + 1]))
            for i in range(len(bounds_v4PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_wavenb_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[0][i * 3 + 2]))
            for i in range(len(bounds_v4PO4[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_wavenb_v4PO4')
        [conffile.write(' {}'.format(bounds_v4PO4[1][i * 3 + 2]))
            for i in range(len(bounds_v4PO4[1]) // 3)]
        conffile.write('\n')
        conffile.write('\n')

        # CO3 area
        conffile.write('CO3bound {0} {1}\n'.format(CO3[0], CO3[1]))
        conffile.write('guess_height_CO3')
        [conffile.write(' {}'.format(guess_CO3[i * 3]))
            for i in range(len(guess_CO3) // 3)]
        conffile.write('\n')
        conffile.write('guess_halffwhm_CO3')
        [conffile.write(' {}'.format(guess_CO3[i * 3 + 1]))
            for i in range(len(guess_CO3) // 3)]
        conffile.write('\n')
        conffile.write('guess_wavenb_CO3')
        [conffile.write(' {}'.format(guess_CO3[i * 3 + 2]))
            for i in range(len(guess_CO3) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_height_CO3')
        [conffile.write(' {}'.format(bounds_CO3[0][i * 3]))
            for i in range(len(bounds_CO3[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_height_CO3')
        [conffile.write(' {}'.format(bounds_CO3[1][i * 3]))
            for i in range(len(bounds_CO3[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_halffwhm_CO3')
        [conffile.write(' {}'.format(bounds_CO3[0][i * 3 + 1]))
            for i in range(len(bounds_CO3[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_halffwhm_CO3')
        [conffile.write(' {}'.format(bounds_CO3[1][i * 3 + 1]))
            for i in range(len(bounds_CO3[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_wavenb_CO3')
        [conffile.write(' {}'.format(bounds_CO3[0][i * 3 + 2]))
            for i in range(len(bounds_CO3[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_wavenb_CO3')
        [conffile.write(' {}'.format(bounds_CO3[1][i * 3 + 2]))
            for i in range(len(bounds_CO3[1]) // 3)]
        conffile.write('\n')
        conffile.write('\n')

        # amide area
        conffile.write('amidebound {0} {1}\n'.format(amide[0], amide[1]))
        conffile.write('guess_height_amide')
        [conffile.write(' {}'.format(guess_amide[i * 3]))
            for i in range(len(guess_amide) // 3)]
        conffile.write('\n')
        conffile.write('guess_halffwhm_amide')
        [conffile.write(' {}'.format(guess_amide[i * 3 + 1]))
            for i in range(len(guess_amide) // 3)]
        conffile.write('\n')
        conffile.write('guess_wavenb_amide')
        [conffile.write(' {}'.format(guess_amide[i * 3 + 2]))
            for i in range(len(guess_amide) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_height_amide')
        [conffile.write(' {}'.format(bounds_amide[0][i * 3]))
            for i in range(len(bounds_amide[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_height_amide')
        [conffile.write(' {}'.format(bounds_amide[1][i * 3]))
            for i in range(len(bounds_amide[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_halffwhm_amide')
        [conffile.write(' {}'.format(bounds_amide[0][i * 3 + 1]))
            for i in range(len(bounds_amide[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_halffwhm_amide')
        [conffile.write(' {}'.format(bounds_amide[1][i * 3 + 1]))
            for i in range(len(bounds_amide[1]) // 3)]
        conffile.write('\n')
        conffile.write('boundsinf_wavenb_amide')
        [conffile.write(' {}'.format(bounds_amide[0][i * 3 + 2]))
            for i in range(len(bounds_amide[0]) // 3)]
        conffile.write('\n')
        conffile.write('boundssup_wavenb_amide')
        [conffile.write(' {}'.format(bounds_amide[1][i * 3 + 2]))
            for i in range(len(bounds_amide[1]) // 3)]
        conffile.write('\n')
