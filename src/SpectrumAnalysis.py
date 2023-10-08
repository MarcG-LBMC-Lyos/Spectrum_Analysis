# LYOS Spectrum Analysis is a program of automatic analysis of bone FTIRM spectrum.
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
#
# You can contact me by mail at : m.gardegaront@gmail.com


# -*- coding: utf-8 -*-
"""
Main File.
Defines the main window (size, Title, Menu, etc.)
"""

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import matplotlib.pyplot as plt
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
import re
import spc


os.chdir(pa.dirname(pa.realpath(__file__)))  # Change working directory to the one containing the main python script.

RESOLUTION = '720x540'  # Original Resolution of the window
SOFTWARE_NAME = 'LYOS Spectrum Analysis'  # Name of the software (appearing in the top bar of the window)
FILENAME = './configure.ini'  # Name of the configuration file for the automatic analysis

matplotlib.use('TkAgg')  # Tells matplotlib to work in a Tk environment


class MainProgram(tk.Tk):

    def __init__(self, *args, **kwargs):
        """ Creates the main Window, the menu bar and launch the main loop """
        tk.Tk.__init__(self, *args, **kwargs) # Tk class initialization

        # Creation of the default configuration file if doesn't exist
        if not os.path.isfile(FILENAME):
            writeDefaultConfig(FILENAME)

        # Definition of the main window's name and geometry
        self.geometry(RESOLUTION)
        self.wm_title(SOFTWARE_NAME)

        # Initialization of attributes
        self.initialize()
        self.readConfigFile()
        self.display = True

        self.function = sf.gaussianSum
        self.blMethod = 'lin'
        self.currentdir = '.'

        self.fig = plt.figure(figsize=(0, 0))
        self.ax = self.fig.add_subplot(111)
        plt.gca().invert_xaxis()
        plt.xlabel('Wave number')
        plt.ylabel('Transmittance / Absorbance')

        # Creation of the plot canvas
        self.plotFrame = tk.Frame(self)
        self.plotFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plotFrame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plotFrame)
        self.toolbar.update()

        # Menu bar creation
        menubar = tk.Menu(self)

        # Menu container
        mainMenus = OrderedDict()

        # File menu
        mainMenus['File'] = OrderedDict()
        mainMenus['File']['Load file (Ctrl-l)'] = self.load
        mainMenus['File']['Save spectrum (Ctrl-s)'] = self.save
        mainMenus['File']['Save config file'] = self.saveConfigFile
        mainMenus['File']['Load config file'] = self.loadConfigFile
        # mainMenus['File']['Configure'] = self.configure
        mainMenus['File']['Quit'] = self._quit

        # Edit menu
        mainMenus['Edit'] = OrderedDict()
        mainMenus['Edit']['Undo (Ctrl-z)'] = self.undo
        mainMenus['Edit']['Redo (Ctrl-y)'] = self.redo
        mainMenus['Edit']['Total spectrum (Undo ROI)'] = self.remROI

        # Preprocessing menu
        mainMenus['Preprocessing'] = OrderedDict()
        mainMenus['Preprocessing']['Spectrum component removal (Ctrl-r)'] = self.spectrumRemoval
        mainMenus['Preprocessing']['Transmittance to Absorbance (Ctrl-a)'] = self.absorbance
        mainMenus['Preprocessing']['Baseline (Ctrl-b)'] = self.baseline
        mainMenus['Preprocessing']['Select ROI (Ctrl-h)'] = self.roiAnalysis
        mainMenus['Preprocessing']['Normalisation to a specific peak (Ctrl-n)'] = self.normalisation

        # Analysis menu
        mainMenus['Analysis'] = OrderedDict()
        mainMenus['Analysis']['Multiple Files Analysis'] = self.massAnalysis
        mainMenus['Analysis']['Automatic Analysis'] = self.autoAnalysis
        mainMenus['Analysis']['Manual Analysis'] = self.manualAnalysis

        # Menu bar organization
        for menu in mainMenus:
            curMenu = tk.Menu(menubar, tearoff=0)
            for subMenu in mainMenus[menu]:
                curMenu.add_command(label=subMenu, command=mainMenus[menu][subMenu])
            menubar.add_cascade(label=menu, menu=curMenu)

        # Display the menu
        self.config(menu=menubar)

        # Binding ctrl-key to the corresponding actions
        self.bind("<Control-z>", self.undo)
        self.bind("<Control-y>", self.redo)
        self.bind("<Control-r>", self.spectrumRemoval)
        self.bind("<Control-a>", self.absorbance)
        self.bind("<Control-b>", self.baseline)
        self.bind("<Control-h>", self.roiAnalysis)
        self.bind("<Control-n>", self.normalisation)
        self.bind("<Control-l>", self.load)
        self.bind("<Control-s>", self.save)

        # Properly quit the application when clicking the cross icon
        self.protocol('WM_DELETE_WINDOW', self._quit)

    def initialize(self):
        self.filepath = None  # Path to the current data file
        self.data_history = []  # History of the changes occured to the datas
        self.data = [None, None]  # Datas extracted from the file (x, y)
        self.curWin = None  # Current window
        self.absorbcount = False  # Determine if absorbance transformation has been done yet
        self.absorbcount_history = []  # History of the changes occured to absorbcount
        self.bounds = [None, None]  # Determine the x boundaries of the view (and only the view)
        self.bounds_history = []  # History of the changes occured to bounds
        self.historyPointer = -1  # Keep track of the position in history when undo action is called
        self.limReached = False
        self.lenLim = 0

    def loadConfigFile(self):
        path_config = askopenfilename(
            title='Load config file',
            initialdir=self.currentdir)
        self.readConfigFile(path=path_config)


    def saveConfigFile(self):
        configParams = {'mmabounds': self.bounds_sub['MMA'],
                        'v1v3PO4': self.areaBounds['v1v3PO4'],
                        'guess_v1v3PO4': self.guess['v1v3PO4'],
                        'bounds_v1v3PO4': self.peakBounds['v1v3PO4'],
                        'v4PO4': self.areaBounds['v4PO4'],
                        'guess_v4PO4': self.guess['v4PO4'],
                        'bounds_v4PO4': self.peakBounds['v4PO4'],
                        'CO3': self.areaBounds['CO3'],
                        'guess_CO3': self.guess['CO3'],
                        'bounds_CO3': self.peakBounds['CO3'],
                        'amide': self.areaBounds['amide'],
                        'guess_amide': self.guess['amide'],
                        'bounds_amide': self.peakBounds['amide']}
        path_config = askopenfilename(
            title='Save config file',
            initialdir=self.currentdir)
        writeDefaultConfig(path_config, writeCurrent=True, **configParams)

    def readConfigFile(self, path=FILENAME):

        self.bounds_sub = None
        self.areaBounds = None
        self.guess = None
        self.peakBounds = None
        self.limit = 10

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
        with open(path, 'r+') as conffile:
            conf = [line.split() for line in conffile]
            for line in conf:
                if not line:
                    continue

                key = 'limit'
                if len(line[0]) >= len(key):
                    if line[0][len(line[0]) - len(key):] == key:
                        self.limit = float(line[1])
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

        self.bounds_sub = bounds_sub
        self.areaBounds = areaBounds
        self.guess = guess
        self.peakBounds = bounds

    def data_change(self, data, **kwargs):
        refresh_ax_only = kwargs.get('refresh_ax_only', False)
        if refresh_ax_only is False:
            self.data = np.array(data)
        if self.bounds[0] is None or self.bounds[1] is None:
            self.bounds = [0, len(self.data[0])]

        self.ax.cla()
        # plt.gca().invert_xaxis()
        self.ax.plot(self.data[0][self.bounds[0]:self.bounds[1]], self.data[1][self.bounds[0]:self.bounds[1]], '-b')
        plt.xlabel('Wave number ($cm^{-1}$)')
        plt.ylabel('Transmittance / Absorbance')
        # recompute the ax.dataLim
        self.ax.relim()
        self.ax.invert_xaxis()
        # update ax.viewLim using the new dataLim
        self.ax.autoscale_view()
        if self.display:
            self.canvas.draw()

        if refresh_ax_only is False:
            self.history()

    def history(self):
        self.historyPointer += 1
        if len(self.data_history) < self.historyPointer+1:
            self.data_history.append(deepcopy(self.data))
            self.absorbcount_history.append(deepcopy(self.absorbcount))
            self.bounds_history.append(deepcopy(self.bounds))
        else:
            self.data_history[self.historyPointer] = deepcopy(self.data)
            self.absorbcount_history[self.historyPointer] = deepcopy(self.absorbcount)
            self.bounds_history[self.historyPointer] = deepcopy(self.bounds)
            for i in range(self.historyPointer+1, len(self.data_history)):
                self.data_history[i] = None
                self.absorbcount_history[i] = None
                self.bounds_history[i] = None

    def undo(self, *args, **kwargs):
        if len(self.data_history) > 1 and self.historyPointer > 0:
            if self.curWin is not None:
                self.curWin._quit()  # Quits current window
                self.curWin = None
            self.historyPointer -= 1
            self.data = self.data_history[self.historyPointer]
            self.absorbcount = self.absorbcount_history[self.historyPointer]
            self.bounds = self.bounds_history[self.historyPointer]

            self.ax.cla()
            self.ax.plot(self.data[0][self.bounds[0]:self.bounds[1]], self.data[1][self.bounds[0]:self.bounds[1]], '-b')
            # recompute the ax.dataLim
            self.ax.relim()
            # update ax.viewLim using the new dataLim
            self.ax.autoscale_view()
            self.canvas.draw()

    def redo(self, *args, **kwargs):
        if len(self.data_history) > self.historyPointer+1:
            if self.data_history[self.historyPointer+1] is None:
                return
            if self.curWin is not None:
                self.curWin._quit()  # Quits current window
                self.curWin = None
            self.historyPointer += 1
            self.data = self.data_history[self.historyPointer]
            self.absorbcount = self.absorbcount_history[self.historyPointer]
            self.bounds = self.bounds_history[self.historyPointer]

            self.ax.cla()
            self.ax.plot(self.data[0][self.bounds[0]:self.bounds[1]], self.data[1][self.bounds[0]:self.bounds[1]], '-b')
            # recompute the ax.dataLim
            self.ax.relim()
            # update ax.viewLim using the new dataLim
            self.ax.autoscale_view()
            self.canvas.draw()

    def load(self, *args, **kwargs):
        """ Load a new .spc file and erase the potential old one """
        if self.curWin is not None:
            self.curWin._quit()  # Quits current window
            self.curWin = None
        self.initialize()  # Reinitialize attributes for new datas

        # Get the file's path or open file explorer pop up if no path is given
        self.filepath = kwargs.get('filepath', None)
        if self.filepath is None:
            try:
                self.filepath = askopenfilename(
                    title='Open .spc spectrum file',
                    filetypes=[('SP and SPC file', '*.sp *.SPC *.spc'), ('SP file', '*.sp *.SP'), ('SPC file', '*.spc'), ('SPC file', '*.SPC'), ('Text file', '*.txt')],
                    initialdir=self.currentdir)
                self.currentdir = os.path.dirname(self.filepath)
                if not self.filepath:
                    self.filepath = None
                    return
            except:
                self.filepath = None
                return

        # Extracting the datas with the spc library
        if self.filepath[-4:] == '.spc' or self.filepath[-4:] == '.SPC':
            f = spc.File(self.filepath)
            x = np.array(f.x)
            y = np.array(f.sub[0].y)
        elif self.filepath[-3:] == '.sp':
            y, x = spload2(self.filepath)
        elif self.filepath[-4:] == '.txt':
            with open(self.filepath, 'r+') as f:
                data = [[float(a.split(' ')[0]), float(a.split(' ')[1])] for a in f]
                x = np.array([data[i][0] for i in range(len(data))])
                y = np.array([data[i][1] for i in range(len(data))])
        xsortarg = np.argsort(x)  # Sorting the x values
        x = x[xsortarg]
        y = y[xsortarg]

        data = [x, y]

        self.arglim = [i for i in range(len(x))]  # Every argument of the points that will be taken in account during the deconvolution

        self.data_change(data)

    def remROI(self, *args, **kwargs):
        firstBound = self.bounds_history[0]
        i = -1
        for bound in self.bounds_history:
            i += 1
            if firstBound != bound:
                break
        for j in range(i, len(self.bounds_history)):
            self.undo()

    def manualAnalysis(self, *args, **kwargs):
        """
        Trigger the manual analysis mode (see ManualAnalysis.py for
        more details)
        """
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creating a ManualAnalysis instance
        self.curWin = ManualAnalysis(self)
        self.curWin.pack(side=tk.TOP)

    def autoAnalysis(self, *args, **kwargs):
        """
        Trigger the automatic analysis mode (see AutomaticAnalysis.py for
        more details)
        """
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creates an AutomaticAnalysis instance
        self.curWin = AutomaticAnalysis(self)
        self.curWin.pack(side=tk.TOP)

    def roiAnalysis(self, *args, **kwargs):
        """
        Trigger the Region of Interest selection mode (see ROIselect.py for
        more details)
        """
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creating a ROIselect instance
        self.curWin = ROIselect(self)
        self.curWin.pack(side=tk.TOP)

    def spectrumRemoval(self, *args, **kwargs):
        """ Triggers the MMA removing algorithm based on the second derivative.
        Use the second derivative to check if the MMA peak is still there. In
        pragmatical words, it minimizes the apparition of inflexion points due
        to the MMA peak.
        Optional **kwargs : "mmafile" (path to the MMA spectrum file) """
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creating a spectrumRem instance
        self.curWin = spectrumRem(self)
        self.curWin.pack(side=tk.TOP)

    def baseline(self, *args, **kwargs):
        """ Find and subtracts the baseline from the datas """
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creating a Baseline instance
        self.curWin = Baseline(self)
        self.curWin.pack(side=tk.TOP)

    def massAnalysis(self, *args, **kwargs):
        """ Automatic analysis of several files """
        self.display = False
        transmcheck = kwargs.get('transmcheck', None)
        mmaremcheck = kwargs.get('mmaremcheck', None)
        baselcheck = kwargs.get('baselcheck', None)
        smoothcheck = kwargs.get('smoothcheck', None)
        filenames = kwargs.get('filenames', None)
        nowarn = kwargs.get('nowarn', None)
        # Asking which parameters to execute
        if transmcheck is None or mmaremcheck is None or baselcheck is None:
            transmcheck = tk.IntVar()
            mmaremcheck = tk.IntVar()
            baselcheck = tk.IntVar()
            smoothcheck = tk.IntVar()
            massParamWin = tk.Toplevel(master=self)
            transm = tk.Checkbutton(master=massParamWin, text="Transmittance to absorbance", variable=transmcheck)
            transm.toggle()
            transm.pack()
            mmarem = tk.Checkbutton(master=massParamWin, text="MMA removal", variable=mmaremcheck)
            mmarem.toggle()
            mmarem.pack()
            basel = tk.Checkbutton(master=massParamWin, text="Baseline", variable=baselcheck)
            basel.toggle()
            basel.pack()
            smooth = tk.Checkbutton(master=massParamWin, text="Smoothing", variable=smoothcheck)
            smooth.toggle()
            smooth.pack()
            okbutt = tk.Button(master=massParamWin, text="OK", command=massParamWin.destroy)
            okbutt.pack()
            massParamWin.bind("<Return>", lambda event: massParamWin.destroy())

            self.wait_window(massParamWin)
            transmcheck = transmcheck.get()
            mmaremcheck = mmaremcheck.get()
            baselcheck = baselcheck.get()
            smoothcheck = smoothcheck.get()
        # Asking the spectra files
        if filenames is None:
            filenames = askopenfilenames(
                title='Open .spc spectrum files',
                filetypes=[('SP file', '*.sp'), ('SPC file', '*.spc'), ('SPC file', '*.SPC'), ('Text file', '*.txt')],
                initialdir=self.currentdir)
        if not filenames:
            print('Error: no file of spectrum given')
            self.display = True
            return
        self.currentdir = os.path.dirname(filenames[0])
        bounds_sub = self.bounds_sub
        if mmaremcheck:
            mmaDatas = {}
            mmafile = kwargs.get('mmafile', None)
            for areaname in bounds_sub:
                if mmafile is None:
                    mmafile = askopenfilename(
                        title='Open .spc {} spectrum file'.format(areaname),
                        filetypes=[('SP file', '*.sp'), ('SPC file', '*.spc'), ('SPC file', '*.SPC')],
                        initialdir=self.currentdir
                        )  # MMA spectrum file
                if not mmafile:
                    print('Warning: no file of spectrum to be subtracted given for the area : ' + areaname)
                    continue
                if mmafile[-4:] == '.spc':
                    mmaF = spc.File(mmafile)
                    mmaDatas[areaname] = [np.array(mmaF.x), np.array(mmaF.sub[0].y)]
                elif mmafile[-3:] == '.sp':
                    y, x = spload2(mmafile)
                    mmaDatas[areaname] = [np.array(x), np.array(y)]
                if transmcheck:
                    mmaDatas[areaname][1] = -np.log10(mmaDatas[areaname][1]/100)

        # Warning window if the Results directory already exists
        path = os.path.dirname(filenames[0])+'/'
        if os.path.exists(os.path.join(path, 'Results')) or os.path.exists(os.path.join(path, 'Results_fail')):
            if not nowarn:
                warn = tk.Toplevel(master=self)
                warn.wm_title('WARNING')
                warningtext = 'WARNING, "Results" directory already exists in:\n' \
                    + path + '. \n\nChange the directory name or move it' + \
                    'somewhere else to avoid erasing its datas'
                warnlab = tk.Label(master=warn, text=warningtext)
                warnlab.pack()
                destroy = tk.Button(master=warn, text='Ok', command=warn.destroy)
                destroy.pack()
                self.wait_window(warn)
            if os.path.exists(os.path.join(path, 'Results')):
                shutil.rmtree(os.path.join(path, 'Results'))
            if os.path.exists(os.path.join(path, 'Results_fail')):
                shutil.rmtree(os.path.join(path, 'Results_fail'))

        # Automatic analysis loop for each files
        """
        # Doesn't work
        waiting = tk.Toplevel(self)
        tk.Label(waiting, text='Wait while your files are being treated...').pack()
        """
        for filename in filenames:
            try:
                print('Analysing file : '+filename)
                # Loads file
                if self.curWin is not None:
                    self.curWin._quit()
                    self.curWin = None
                self.load(filepath=filename)
                if filename is None:
                    continue

                # Transmittance to absorbance
                if transmcheck:
                    self.absorbance()

                # Check if absorbance threshold is exceeding
                self.arglim = []
                self.lenLim = 0
                self.limReached = False
                argnotlim = []
                for i in range(len(self.data[1])):
                    if self.data[0][i] <= self.areaBounds['v1v3PO4'][1] and self.data[0][i] >= self.areaBounds['v1v3PO4'][0]:
                        if self.data[1][i] <= self.limit:
                            self.arglim.append(i)
                        else:
                            self.limReached = True
                            argnotlim.append(i)

                # MMA removal
                if mmaremcheck:
                    for areaname in mmaDatas:
                        self.data_change(removal(self.data, mmaDatas[areaname],
                                         bounds_sub[areaname]))

                # Baseline subtraction
                if baselcheck:
                    bl = sf.baseline(self.data[0], self.data[1], method=self.blMethod)
                    self.data_change([self.data[0], self.data[1] - bl])

                # Upscaling and correcting the curve if it reaches the limit
                if smoothcheck:
                    xold = deepcopy(self.data[0])
                    yold = deepcopy(self.data[1])

                    try:
                        self.lenLim = abs(self.data[0][argnotlim[0]] - self.data[0][argnotlim[-1]])
                    except:
                        self.lenLim = 0
                        pass

                    finterp = interp1d(self.data[0][self.arglim], self.data[1][self.arglim], kind='cubic')
                    x1 = abs(self.data[0] - self.areaBounds['v1v3PO4'][0]).argmin()
                    x2 = abs(self.data[0] - self.areaBounds['v1v3PO4'][1]).argmin()
                    x = np.linspace(self.data[0][x1], self.data[0][x2], 5000)
                    y = finterp(x)
                    x = list(self.data[0][0:x1]) + list(x) + list(self.data[0][x2+1:])
                    y = list(self.data[1][0:x1]) + list(y) + list(self.data[1][x2+1:])

                    self.bounds = [0, len(x)]
                    self.data_change([x, y])

                    self.ax.plot(xold[self.arglim], yold[self.arglim], '--r')
                    self.ax.plot(xold, yold, '-g')

                # Automatic analysis
                self.autoAnalysis()
                self.curWin.refreshCallback()  # Deconvolution
                self.curWin.saveCallback()  # Calculating and saving the results
                print('Done.')
            except Exception as e:
                print('ERROR FILE : {}'.format(filename))
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)

        """waiting.destroy()"""
        #########################################
        # GLOBAL RESULT FILE
        dirPaths = [os.path.join(path, 'Results')]
        mean = [[] for i in range(8)]
        std = [[] for i in range(8)]
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

            # Sort function
            def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
                return [int(text) if text.isdigit() else text.lower()
                        for text in re.split(_nsre, s)]
            dataFiles = sorted(dataFiles, key=natural_sort_key)

            results = [[] for i in range(8)]
            for fname in dataFiles:
                with open(pa.join(dirPath, fname), 'r+') as f:
                    f.readline()
                    vals = f.readline().split(',')
                    if vals[0] == '\n':
                        vals = f.readline().split(',')
                    print(vals)
                    results[0].append(float(vals[0]))
                    results[1].append(float(vals[1]))
                    results[2].append(float(vals[2]))
                    results[3].append(float(vals[3]))
                    results[4].append(float(vals[4]))
                    results[5].append(float(vals[5]))
                    results[6].append(float(vals[6]))
                    results[7].append(float(vals[7]))
                    name2.append(fname)

            mean = results

        print(dataFiles)

        meantemp = [[] for i in range(8)]
        for j in range(len(name2)):
            for i in range(len(meantemp)):
                meantemp[i].append(mean[i][j])

        mean = meantemp

        arg = np.argsort(mean[0])

        mean = np.array(mean)

        label = ['Filename', 'Mineralization Index', 'Carbonation',
                 'Crystallinity', 'Mineral Maturity', 'Collagen Maturity', 'Collagen Maturity (intensity)', 'Mineralization Index (intensity)', 'Mineral Maturity (intensity)', 'v1v3PO4 area', 'Cut length']
        fileSave = Workbook(os.path.join(path, 'Results')+'/Results.xlsx')
        worksheet = fileSave.add_worksheet()

        worksheet.write_row(0, 0, label)
        for i in range(len(dataFiles)):
            try:
                worksheet.write_row(i + 1, 0, [dataFiles[i], mean[0][i], mean[1][i], mean[2][i], mean[3][i], mean[4][i],
                                               mean[5][i], mean[6][i], mean[7][i]])
            except:
                worksheet.write_row(i + 1, 0,
                                    [dataFiles[i], str(mean[0][i]), str(mean[1][i]), str(mean[2][i]), str(mean[3][i]),
                                     str(mean[4][i]), str(mean[5][i]), str(mean[6][i]), str(mean[7][i])])

        worksheet.write_row(len(dataFiles)+2, 0, ['AVERAGE','=AVERAGE(B2:B'+str(len(dataFiles)+1)+')','=AVERAGE(C2:C'+str(len(dataFiles)+1)+')','=AVERAGE(D2:D'+str(len(dataFiles)+1)+')','=AVERAGE(E2:E'+str(len(dataFiles)+1)+')','=AVERAGE(F2:F'+str(len(dataFiles)+1)+')','=AVERAGE(G2:G'+str(len(dataFiles)+1)+')','=AVERAGE(H2:H'+str(len(dataFiles)+1)+')','=AVERAGE(I2:I'+str(len(dataFiles)+1)+')'])



        fileSave.close()

        # SAVING USED CONFIG FILE
        try:
            configParams = {'mmabounds': self.bounds_sub['MMA'],
                            'v1v3PO4': self.areaBounds['v1v3PO4'],
                            'guess_v1v3PO4': self.guess['v1v3PO4'],
                            'bounds_v1v3PO4': self.peakBounds['v1v3PO4'],
                            'v4PO4': self.areaBounds['v4PO4'],
                            'guess_v4PO4': self.guess['v4PO4'],
                            'bounds_v4PO4': self.peakBounds['v4PO4'],
                            'CO3': self.areaBounds['CO3'],
                            'guess_CO3': self.guess['CO3'],
                            'bounds_CO3': self.peakBounds['CO3'],
                            'amide': self.areaBounds['amide'],
                            'guess_amide': self.guess['amide'],
                            'bounds_amide': self.peakBounds['amide']}
            writeDefaultConfig(os.path.join(path, 'Results')+'/Configure.ini', writeCurrent=True, **configParams)
        except:
            print("\nWARNING: COULD NOT SAVE THE CONFIG FILE!\n")

        ###########################################

        self.display = True

    def absorbance(self, *args, **kwargs):
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        if self.absorbcount is True:
            return
        self.absorbcount = True
        data = [self.data[0], -np.log10(self.data[1]/100)]
        self.data_change(data)

        self.mmaabsorb_check = True

    def smoothing(self, *args, **kwargs):
        nbpoints = kwargs.get('nbpoints', 10000)
        finterp = interp1d(self.data[0][self.arglim], self.data[1][self.arglim], kind='cubic')
        x = np.linspace(min(self.data[0]), max(self.data[0]), nbpoints)
        y = finterp(x)
        self.bounds = [0, len(x)]
        self.data_change([x, y])

    def normalisation(self, *args, **kwargs):
        if self.filepath is None:  # If no file loaded
            self.load()
            # If failed to load, abort the execution
            if self.filepath is None:
                return

        if self.curWin is not None:
            self.curWin._quit()
            self.curWin = None

        # Creating a Baseline instance
        self.curWin = Normalisation(self)
        self.curWin.pack(side=tk.TOP)

    def configure(self, *args, **kwargs):
        Configure(self)

    def save(self, *args, **kwargs):
        global name
        f = self.f
        f.x = self.x
        f.sub[0].y = self.y
        direc = askdirectory(initialdir=self.currentdir)
        name = 'Spectrum'
        nametop = tk.Toplevel(master=self)
        tk.Label(master=nametop, text='File''s name').pack(side=tk.LEFT)
        name_edit = tk.Entry(master=nametop)
        name_edit.pack(side=tk.LEFT)

        def save_callback():
            global name
            name = name_edit.get()
            nametop.destroy()

        ok_butt = tk.Button(master=nametop, text='Save', command=save_callback)
        ok_butt.pack(side=tk.LEFT)

        self.wait_window(nametop)
        f.write_file(direc+'/'+name+'.spc')

    def _quit(self, *args, **kwargs):
        """ Properly quit the program """
        self.quit()
        self.destroy()


"""
# Begin the program
dirnames = ['C:/Users/Delphine/Desktop/Marc/MULTIPS/Lateral/'+str(i) for i in range(1, 29) if i != 7 and i != 9 and i != 10 and i != 16 and i != 19 and i != 22]
# ['C:/Users/Delphine/Desktop/Marc/MULTIPS/Medial/'+str(i) for i in range(33, 58)]
"""
['C:/Users/Delphine/Desktop/Marc/Temoin_crete_illiac/025_001_120M0',
            'C:/Users/Delphine/Desktop/Marc/Temoin_crete_illiac/025_001_130M0',
            'C:/Users/Delphine/Desktop/Marc/Temoin_crete_illiac/025_003_1115M0',
            'C:/Users/Delphine/Desktop/Marc/Temoin_crete_illiac/025_003_1117M0']
"""

program = MainProgram()
program.update_idletasks()
program.update()

for folder in dirnames:
    if os.path.exists(os.path.join(folder, 'Results')):
        shutil.rmtree(os.path.join(folder, 'Results'))
    filenamestemp = [filename for filename in os.listdir(folder) if filename.endswith('.SPC')]
    filenames = []
    mmafile = None
    for filename in filenamestemp:
        if filename[:3] == 'MMA':
            mmafile = folder+'/'+filename
        else:
            filenames.append(folder+'/'+filename)

    makeys = {'transmcheck':True,
              'mmaremcheck':True,
              'baselcheck':True,
              'smoothcheck':True,
              'filenames':filenames,
              'mmafile':mmafile,
              'nowarn':True}
    program.massAnalysis(**makeys)

program._quit()

#"""

#"""# Begin the program
program = MainProgram()
tk.mainloop()
#"""
