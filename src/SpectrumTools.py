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


# -*- coding: utf-8 -*-

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2TkAgg
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import matplotlib.pyplot as plt
import spc
import spectrumFit as sf
from ManualAnalysis import ManualAnalysis
from AutomaticAnalysis import AutomaticAnalysis
from ROIselect import ROIselect
import tkinter as tk
import numpy as np
import os
from Configure import Configure
from writeDefaultConfig import writeDefaultConfig
from collections import OrderedDict
from copy import deepcopy
from spectrumRem import spectrumRem, removal
from Baseline import Baseline
from Normalisation import Normalisation
import sys
import shutil
from scipy.interpolate import interp1d
from spload import spload2
import os.path as pa
from xlsxwriter import Workbook
import csv
from threading import Thread


os.chdir(pa.dirname(pa.realpath(__file__)))  # Change working directory to the one containing the main python script.

def load(filepath):
    """ Load a new .spc file and erase the potential old one """
    # Extracting the datas with the spc library
    if filepath[-4:] == '.spc' or filepath[-4:] == '.SPC':
        f = spc.File(filepath)
        x = np.array(f.x)
        y = np.array(f.sub[0].y)
    elif filepath[-3:] == '.sp':
        y, x = spload2(filepath)
    elif filepath[-4:] == '.txt':
        with open(filepath, 'r+') as f:
            data = [[float(a.split(' ')[0]), float(a.split(' ')[1])] for a in f]
            x = np.array([data[i][0] for i in range(len(data))])
            y = np.array([data[i][1] for i in range(len(data))])
    xsortarg = np.argsort(x)  # Sorting the x values
    x = x[xsortarg]
    y = y[xsortarg]

    data = np.array([x, y])

    return data


def readConfigFile(configFile):
    bounds_sub = None
    areaBounds = None
    guess = None
    peakBounds = None
    limit = 10

    areaname = []  # Areas' name
    bounds_sub = {}  # Boundaries for each peak used to subtract a spectrum (e.g. MMA) (by name)
    areaBounds = {}  # Areas' boundaries (by name)
    guess_height = {}  # Initial guess for each peak's height (by area's name)
    guess_halffwhm = {}  # Initial guess for each peak's half FWHM (by area's name)
    guess_wavenb = {}  # Initial guess for each peak's wave number (by area's name)
    boundsinf_height = {}  # Inferior boundary for each peak's height (by area's name)
    boundsinf_halffwhm = {}  # Inferior boundary for each peak's half FWHM (by area's name)
    boundsinf_wavenb = {}  # Inferior boundary for each peak's wave number (by area's name)
    boundssup_height = {}  # Superior boundary for each peak's height (by area's name)
    boundssup_halffwhm = {}  # Superior boundary for each peak's half FWHM (by area's name)
    boundssup_wavenb = {}  # Superior boundary for each peak's wave number (by area's name)

    # Reading the configurations' file
    with open(configFile, 'r+') as conffile:
        conf = [line.split() for line in conffile]
        for line in conf:
            if not line:
                continue

            key = 'limit'
            if len(line[0]) >= len(key):
                if line[0][len(line[0]) - len(key):] == key:
                    limit = float(line[1])
                    continue

            key = 'bound_sub'
            if len(line[0]) >= len(key):
                if line[0][len(line[0]) - len(key):] == key:
                    bounds_sub[line[0][0:-len(key)]] = [float(val) for val in line[1:]]
                    continue

            key = 'bound'
            if len(line[0]) >= len(key):
                if line[0][len(line[0]) - len(key):] == key:
                    areaname.append(line[0][0:-len(key)])
                    areaBounds[line[0][0:-len(key)]] = [float(val) for val in line[1:]]
                    continue

            key = 'guess_height'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    guess_height[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'guess_halffwhm'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    guess_halffwhm[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'guess_wavenb'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    guess_wavenb[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundsinf_height'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundsinf_height[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundsinf_halffwhm'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundsinf_halffwhm[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundsinf_wavenb'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundsinf_wavenb[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundssup_height'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundssup_height[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundssup_halffwhm'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundssup_halffwhm[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

            key = 'boundssup_wavenb'
            if len(line[0]) >= len(key):
                if line[0][0:len(key)] == key:
                    boundssup_wavenb[line[0][len(key)+1:]] = [float(val) for val in line[1:]]
                    continue

        guess = {}
        bounds = {}
        for j in areaname:
            guess[j] = np.zeros(len(guess_height[j])*3)
            boundsinf = np.zeros(len(guess_height[j])*3)
            boundssup = np.zeros(len(guess_height[j])*3)
            for i in range(len(guess_height[j])):
                guess[j][i*3] = guess_height[j][i]
                guess[j][i*3+1] = guess_halffwhm[j][i]
                guess[j][i*3+2] = guess_wavenb[j][i]
                boundsinf[i*3] = boundsinf_height[j][i]
                boundsinf[i*3+1] = boundsinf_halffwhm[j][i]
                boundsinf[i*3+2] = boundsinf_wavenb[j][i]
                boundssup[i*3] = boundssup_height[j][i]
                boundssup[i*3+1] = boundssup_halffwhm[j][i]
                boundssup[i*3+2] = boundssup_wavenb[j][i]
            bounds[j] = [boundsinf, boundssup]

    return [bounds_sub, areaBounds, guess, bounds]


def deconvolute(filepath, data, areaBounds, areaBoundsId, blMethod, guess,
                    peakBounds, limReached, lenLim, xLim, data_orig):
    """ Automatic deconvolution callback """
    poptsave = {}
    for area in areaBounds:  # For every area
        # Restraining the ROI to the analysed area
        # Getting Id of the boundaries
        x = data[0][areaBoundsId[area][0]:areaBoundsId[area][1]+1]
        y = data[1][areaBoundsId[area][0]:areaBoundsId[area][1]+1]

        # Calculating and subtracting the baseline of the ROI
        bl = sf.baseline(x, y, method=blMethod)
        y = y - bl

        # Bone specific
        if area == 'amide':
            # Height constraint of the near 1600 peak
            try:
                heights = np.array([peakBounds['amide'][1][i*3] for i in range(int(len(peakBounds['amide'][1])/3))])
                wavenbs = np.array([peakBounds['amide'][1][i*3+2] for i in range(int(len(peakBounds['amide'][1])/3))])
                arg1600 = abs(wavenbs - 1600).argmin()
                arg1600x = abs(x - 1600).argmin()
                peakBounds['amide'][1][arg1600*3] = y[arg1600x]*2./3
            except Exception as e:
                print(str(e))

        if area == 'v1v3PO4':
            # Height constraint of the near 1060 peak (based on the 1030 peak)
            try:
                heights = np.array([peakBounds['v1v3PO4'][1][i*3] for i in range(int(len(peakBounds['v1v3PO4'][1])/3))])
                wavenbs = np.array([peakBounds['v1v3PO4'][1][i*3+2] for i in range(int(len(peakBounds['v1v3PO4'][1])/3))])
                arg1030x = abs(x - 1030).argmin()
                arg1060 = abs(wavenbs - 1060).argmin()
                peakBounds['v1v3PO4'][1][arg1060*3] = y[arg1030x]*1./3

                arg1030 = abs(wavenbs - 1030).argmin()
                arg1110 = abs(wavenbs - 1110).argmin()
                arg1110x = abs(x - 1110).argmin()
                peakBounds['v1v3PO4'][0][arg1030*3] = y[arg1030x]*4./5
                peakBounds['v1v3PO4'][0][arg1110*3] = y[arg1110x]*4./5
                guess['v1v3PO4'][arg1030*3] = y[arg1030x]*4./5
                guess['v1v3PO4'][arg1110*3] = y[arg1110x]*4./5

            except Exception as e:
                print(str(e))

        popt = sf.spectrumFit(sf.gaussianSum,
                              x, y, guess=guess[area],
                              bounds=peakBounds[area])

        # Saving the ROI peaks
        poptsave[area] = popt


    # v1v3PO4 parameters
    for area in areaBounds:
        if area == 'v1v3PO4':
            # Restraining datas to the ROI
            x = data[0][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            y = data[1][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            bl = sf.baseline(x, y, method=blMethod)
            y = y - bl

            # Calculating the area/integral of the concerned peaks
            # Resampling the datas for better area calculation

            v1v3PO4_area = abs(np.trapz(y, x))
            if limReached is True:
                v1v3PO4_area *= (1 - np.exp(0.03878*lenLim - 6.652))

            wavenbs = np.array([poptsave[area][j * 3 + 2] for j in range(int(len(poptsave[area]) / 3))])

            peak1030id_intensity = abs(x - 1030).argmin()
            peak1110id_intensity = abs(x - 1110).argmin()

            peak1030_intensity = y[peak1030id_intensity]
            peak1110_intensity = y[peak1110id_intensity]

            peak1030id = abs(wavenbs - 1030).argmin()
            peak1110id = abs(wavenbs - 1110).argmin()
            x = np.linspace(-500, 500, 10000)
            peak1030_area = abs(
                np.trapz(
                    sf.gaussianSum(
                        x, poptsave[area][peak1030id * 3], poptsave[area][peak1030id * 3 + 1], 0),
                    x))
            peak1110_area = abs(
                np.trapz(
                    sf.gaussianSum(
                        x, poptsave[area][peak1110id * 3], poptsave[area][peak1110id * 3 + 1], 0),
                    x))


            #print peak1030_area
            #print self.poptsave[i][peak1030id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1030id * 3+1]/np.sqrt(np.log(2)))**2)
            #print peak1110_area
            #print self.poptsave[i][peak1110id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1110id * 3+1]/np.sqrt(np.log(2)))**2)

    # v4PO4 parameters
        if area == 'v4PO4':
            # Restraining datas to the ROI
            x = data[0][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            y = data[1][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            bl = sf.baseline(x, y, method=blMethod)
            y = y - bl

            # Calculating the fwhm of the concerned peaks
            wavenbs = np.array([poptsave[area][j * 3 + 2] for j in range(int(len(poptsave[area]) / 3))])

            peak604id = abs(wavenbs - 604).argmin()
            xtest = np.linspace(0, 100, 2000)
            peak604 = np.array(sf.gaussianSum(xtest, poptsave[area][peak604id * 3],
                                              poptsave[area][peak604id * 3 + 1], 0))
            peak604_fwhm = xtest[abs(peak604 - max(peak604)/2).argmin()] * 2

    # CO3 parameters
        if area == 'CO3':
            # Restraining datas to the ROI
            x = data[0][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            y = data[1][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            yorig = y
            xorig = x

            idmin = len(y) - 1
            y = yorig[:idmin+1]
            x = xorig[:idmin+1]
            idstart = 0
            a = (yorig[idmin] - yorig[idstart])/(xorig[idmin] - xorig[idstart])
            b = yorig[0] - a*xorig[0]
            y = y - (a*x+b)
            xstand = x
            ystand = y
            while min(y) < 0:
                idmin -= 1
                y = yorig[:idmin+1]
                x = xorig[:idmin+1]
                a = (yorig[idmin] - yorig[idstart])/(xorig[idmin] - xorig[idstart])
                b = yorig[0] - a*xorig[0]
                y = y - (a*x+b)
                if x[-1] < 887:
                    x = xstand
                    y = ystand
                    break

            # Calculating the area/integral of the concerned peaks
            CO3_area = np.trapz(y, x)

    # amide parameters
        if area == 'amide':
            # Restraining datas to the ROI
            x = data[0][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            y = data[1][areaBoundsId[area][0]:
                                    areaBoundsId[area][1]+1]
            bl = sf.baseline(x, y, method=blMethod)
            y = y - bl

            # Calculating the area/integral of the concerned peaks
            wavenbs = np.array([poptsave[area][j * 3 + 2] for j in range(int(len(poptsave[area]) / 3))])

            peak1660id = abs(wavenbs - 1660).argmin()
            peak1690id = abs(wavenbs - 1690).argmin()
            x = np.linspace(-500, 500, 10000)
            peak1660_area = abs(
                np.trapz(
                    sf.gaussianSum(
                        x, poptsave[area][peak1660id * 3], poptsave[area][peak1660id * 3 + 1], 0),
                    x))
            peak1690_area = abs(
                np.trapz(
                    sf.gaussianSum(
                        x, poptsave[area][peak1690id * 3], poptsave[area][peak1690id * 3 + 1], 0),
                    x))
            #print peak1690_area
            #print self.poptsave[i][peak1690id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1690id * 3+1]/np.sqrt(np.log(2)))**2)
            #print peak1660_area
            #print self.poptsave[i][peak1660id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1660id * 3+1]/np.sqrt(np.log(2)))**2)

    # Amide I area
    try:
        x1 = 1592
        x2 = 1730
        x1id = abs(data[0] - x1).argmin()
        x2id = abs(data[0] - x2).argmin()
        x = data[0][min(x1id, x2id):max(x1id, x2id)+1]
        y = data[1][min(x1id, x2id):max(x1id, x2id)+1]
        a = (y[-1] - y[0])/(x[-1] - x[0])
        b = y[-1] - a*x[-1]
        bl = a*x+b

        y = y - bl

        amideI_area = np.trapz(y, x)
    except:
        pass

    # Caracterisation parameters calculation
    failure = False
    try:
        MI = v1v3PO4_area / amideI_area
        if MI < 1.96 or MI > 8.56:
            failure = True
    except:
        MI = None
        pass
    try:
        carbonation = CO3_area / v1v3PO4_area
        if carbonation < 0.00301 or carbonation > 0.0120:
            failure = True
    except:
        carbonation = None
        pass
    try:
        crystallinity = peak604_fwhm
        if crystallinity < 23.31 or crystallinity > 29.81:
            failure = True
    except:
        crystallinity = None
        pass
    try:
        MM = peak1030_area / peak1110_area
        MM_intensity = peak1030_intensity / peak1110_intensity
        if MM < 0.145 or MM > 3.675:
            failure = True
    except:
        MM = None
        MM_intensity = None
        pass
    try:
        CM = peak1660_area / peak1690_area
        if CM < 1.83 or CM > 6.58:
            failure = True
    except:
        CM = None
        pass

    # Computing the 1030/1110 ratio calculated on the global spectrum
    """
    try:
        area = 'v1v3PO4'
        matMin = deepcopy(self.poptsave[area])

    except:
    """
    # Exporting the results as a csv file
    # First line of the csv file
    label = ['Mineralization Index', 'Carbonation',
             'Crystallinity', 'Mineral Maturity', 'Collagen Maturity', 'Mineral Maturity (intensity)', 'v1v3PO4 area', 'Cut length']
    # Second line
    results = [MI, carbonation, crystallinity, MM, CM, MM_intensity, v1v3PO4_area, lenLim]
    # Extracting the original filname
    path = os.path.dirname(filepath)
    name = os.path.splitext(os.path.basename(filepath))[0]
    # Output file name
    filename = path + '/Results/' + name + '_results.csv'
    filenameFail = path + '/Results_fail/' + name + '_results.csv'
    if not os.path.exists(path + '/Results/'):
        os.makedirs(path + '/Results/')
    #if failure and not os.path.exists(path + '/Results_fail/'):
        #os.makedirs(path + '/Results_fail/')
    # Writing the filename
    #f not failure:

    with open(filename, 'w') as myfile:
        c = csv.writer(myfile)
        c.writerow(label)
        c.writerow(results)
    #else:
        #with open(filenameFail, 'w') as myfile:
            #c = csv.writer(myfile)
            #c.writerow(label)
            #c.writerow(results)

    # Saving the plot
    xmin = min([areaBounds[area][0] for area in areaBounds])
    xmax = max([areaBounds[area][1] for area in areaBounds])
    totRange = xmax-xmin
    xmin -= 0.02*totRange
    xmax += 0.02*totRange

    colorArea = ['r', 'g', 'm', 'c', 'y']
    colorPeak = ['#1dba27', '#dc9100', '#9d27d2', '#1a5bb4', '#ddd526']
    n = 0
    plt.figure()
    plt.xlim(xmin, xmax)
    plt.gca().invert_xaxis()
    plt.plot(data_orig[0], data_orig[1], '-b')
    for area in areaBounds:
        x = data[0][areaBoundsId[area][0]:
                                areaBoundsId[area][1]+1]
        y = data[1][areaBoundsId[area][0]:
                                areaBoundsId[area][1]+1]
        bl = sf.baseline(x, y, method=blMethod)
        for i in range(int(len(poptsave[area]/3))):
            plt.plot(x, sf.gaussianSum(x, *poptsave[area][i*3:i*3+3])+bl, colorPeak[i%len(colorPeak)])
        if area == 'v1v3PO4' and xLim is not None:
            argxlim = [abs(x - xLim[0]).argmin(), abs(x - xLim[1]).argmin()]
            plt.plot(x[argxlim], y[argxlim], 'r--')
            print(x[argxlim])
        plt.plot(x, sf.gaussianSum(x, *poptsave[area])+bl, 'k')
        plt.plot(x, bl, colorArea[n%len(colorArea)])
        plt.axvline(x[0], color=colorArea[n%len(colorArea)], linestyle='--', label=area)
        plt.axvline(x[-1], color=colorArea[n%len(colorArea)], linestyle='--')
        n += 1
    plt.legend()

    origsize = plt.gcf().get_size_inches()
    plt.gcf().set_size_inches(origsize[0]*3, origsize[1])
    #if not failure:
    plt.savefig(path + '/Results/' + name + '_results_spectrum.png', dpi=300)
    #else:
        #plt.savefig(path + '/Results_fail/' + name + '_results_spectrum.png', dpi=300)
    plt.gcf().set_size_inches(origsize[0], origsize[1])

    return poptsave


def absorbance(data):
    data = np.array([data[0], -np.log10(data[1]/100)])
    return data


def massAnalysis(filenames, configFile, mmafile=None, bounds_sub=[1706, 1826],
                 transmcheck=True, mmaremcheck=True, baselcheck=True,
                 smoothcheck=True, limit=1.5):

    bounds_sub, areaBounds, guess, bounds = readConfigFile(configFile)
    bounds_sub = list(bounds_sub.values())[0]
    if mmafile is None:
        mmaremcheck = False

    if mmaremcheck:
        mmaDatas = []
        if mmafile[-4:] == '.spc':
            mmaF = spc.File(mmafile)
            mmaDatas = [np.array(mmaF.x), np.array(mmaF.sub[0].y)]
        elif mmafile[-3:] == '.sp':
            y, x = spload2(mmafile)
            mmaDatas = [np.array(x), np.array(y)]
        if transmcheck:
            mmaDatas = absorbance(mmaDatas)

    # Warning window if the Results directory already exists
    path = os.path.dirname(filenames[0])+'/'
    if os.path.exists(os.path.join(path, 'Results')):
        shutil.rmtree(os.path.join(path, 'Results'))

    # Automatic analysis loop for each files
    for filename in filenames:
        try:
            print('Analysing file : '+filename)
            # Loads file
            data = np.array(load(filename))

            # Transmittance to absorbance
            if transmcheck:
                data = np.array(absorbance(data))

            # Check if absorbance threshold is exceeding
            arglim = []
            lenLim = 0
            limReached = False
            argnotlim = []
            for i in range(len(data[1])):
                if data[0][i] <= areaBounds['v1v3PO4'][1] and data[0][i] >= areaBounds['v1v3PO4'][0]:
                    if data[1][i] <= limit:
                        arglim.append(i)
                    else:
                        limReached = True
                        argnotlim.append(i)

            if argnotlim:
                xlim = [data[0][argnotlim[0]], data[0][argnotlim[-1]]]
            else:
                xlim = None

            # MMA removal
            if mmaremcheck:
                data = np.array(removal(data, mmaDatas, bounds_sub))

            # Baseline subtraction
            if baselcheck:
                bl = sf.baseline(data[0], data[1], method='lin')
                data = np.array([data[0], data[1] - bl])

            data_orig = deepcopy(data)
            # Upscaling and correcting the curve if it reaches the limit
            if smoothcheck:
                xold = deepcopy(data[0])
                yold = deepcopy(data[1])

                try:
                    lenLim = abs(data[0][argnotlim[0]] - data[0][argnotlim[-1]])
                except:
                    lenLim = 0
                    pass

                finterp = interp1d(data[0][arglim], data[1][arglim], kind='cubic')
                x1 = abs(data[0] - areaBounds['v1v3PO4'][0]).argmin()
                x2 = abs(data[0] - areaBounds['v1v3PO4'][1]).argmin()
                print(x1, x2)
                print(arglim[0], arglim[-1])
                x = np.linspace(data[0][x1], data[0][x2], 5000)
                y = finterp(x)
                x = list(data[0][0:x1]) + list(x) + list(data[0][x2+1:])
                y = list(data[1][0:x1]) + list(y) + list(data[1][x2+1:])

                data = np.array([x, y])

            # Automatic analysis
            areaBoundsId = {}
            for area in areaBounds:
                x1id = abs(data[0] - areaBounds[area][0]).argmin()
                x2id = abs(data[0] - areaBounds[area][1]).argmin()
                areaBoundsId[area] = sorted([x1id, x2id])

            deconvolute(filename, data, areaBounds, areaBoundsId, 'lin', guess, bounds, limReached, lenLim, xlim, data_orig)
            print('Done.')
        except:
            print('ERROR DECONVOLUTION : file %s'%(filename))

    #########################################
    # GLOBAL RESULT FILE
    dirPaths = [os.path.join(path, 'Results')]
    mean = [[] for i in range(6)]
    std = [[] for i in range(6)]
    i=0
    name2 = []
    for dirPath in dirPaths:
        dataFilesTemp = [f for f in os.listdir(dirPath) if pa.isfile(pa.join(dirPath, f))]
        dataFiles = []
        for f in dataFilesTemp:
            crit = '.csv'
            if len(f) > len(crit):
                if f[-len(crit):] == '.csv':
                    dataFiles.append(f)
        #dataFiles.sort(key=sortfunc)

        results = [[] for i in range(6)]
        for fname in dataFiles:
            with open(pa.join(dirPath, fname), 'r+') as f:
                f.readline()
                vals = f.readline().split(',')
                if vals[0] == '\n':
                    vals = f.readline().split(',')
                results[0].append(float(vals[0]))
                results[1].append(float(vals[1]))
                results[2].append(float(vals[2]))
                results[3].append(float(vals[3]))
                results[4].append(float(vals[4]))
                results[5].append(float(vals[5]))
                name2.append(fname)

        for i in range(5):
            mean = results

    print(dataFiles)

    meantemp = [[], [], [], [], [], []]
    for j in range(len(name2)):
        meantemp[0].append(mean[0][j])
        meantemp[1].append(mean[1][j])
        meantemp[2].append(mean[2][j])
        meantemp[3].append(mean[3][j])
        meantemp[4].append(mean[4][j])
        meantemp[5].append(mean[5][j])

    mean = meantemp

    arg = np.argsort(mean[0])

    mean = np.array(mean)

    label = ['Mineralization Index', 'Carbonation',
             'Crystallinity', 'Mineral Maturity', 'Collagen Maturity', 'Mineral Maturity (intensity)']
    fileSave = Workbook(os.path.join(path, 'Results')+'/Results.xlsx')
    worksheet = fileSave.add_worksheet()

    worksheet.write_row(0, 0, label)
    for i in range(len(dataFiles)):
        worksheet.write_row(i+1, 0, [dataFiles[i],str(mean[0][i]),str(mean[1][i]),str(mean[2][i]),str(mean[3][i]),str(mean[4][i]),str(mean[5][i])])
    worksheet.write_row(len(dataFiles)+2, 0, ['AVERAGE',str(np.mean(mean[0])),str(np.mean(mean[1])),str(np.mean(mean[2])),str(np.mean(mean[3])),str(np.mean(mean[4])),str(np.mean(mean[5]))])

    fileSave.close()


filenames = askopenfilenames(
    title='Open .spc spectrum files',
    filetypes=[('SP file', '*.sp'), ('SPC file', '*.spc'), ('SPC file', '*.SPC'), ('Text file', '*.txt')])
mmafile = askopenfilename(
    title='Open .spc MMA spectrum file',
    filetypes=[('SP file', '*.sp'), ('SPC file', '*.spc'), ('SPC file', '*.SPC'), ('Text file', '*.txt')])
configFile = 'D:/Marc2/MULTIPS/FTIRM/Spectrum_Analysis_new/configure.ini'

massAnalysis(filenames, configFile, mmafile=mmafile, bounds_sub=[1706, 1826],
                 transmcheck=True, mmaremcheck=True, baselcheck=True,
                 smoothcheck=True, limit=1.5)



"""
class MassAn(Thread):

    def __init__(self, filenames, mmafile, configFile):
        Thread.__init__(self)
        self.filenames = filenames
        print(len(filenames))
        self.mmafile = mmafile
        self.configFile = configFile

    def run(self):
        massAnalysis(filenames, configFile, mmafile=mmafile, bounds_sub=[1706, 1826],
                         transmcheck=True, mmaremcheck=True, baselcheck=True,
                         smoothcheck=True, limit=1.5)


filenames = ['D:/Marc2/MULTIPS/FTIRM/Lateral_Thinner/1_precision/map000'+str(i)+'.sp' for i in range(100, 301)]
mmafile = 'D:/Marc2/MULTIPS/FTIRM/Lateral_Thinner/1_precision/MMA1.sp'
configFile = 'D:/Marc2/MULTIPS/FTIRM/Spectrum_Analysis_new/configure.ini'

nbCore = 4

threads = []
for i in range(nbCore):
    start = int(i*len(filenames)/nbCore)
    end = int((i+1)*len(filenames)/nbCore)
    threads.append(MassAn(filenames[start:end], mmafile, configFile))

for i in range(nbCore):
    threads[i].start()

for i in range(nbCore):
    threads[i].join()
"""
