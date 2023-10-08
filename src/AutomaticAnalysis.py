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

"""
Automatic Analysis class.
Defines an automatic analysis window and tools.
"""

import numpy as np
import matplotlib.pyplot as plt
import spectrumFit as sf
import tkinter as tk
from tkinter.ttk import Notebook
from functools import partial
import csv
import os.path
from copy import deepcopy
from random import uniform

PICKERSIZE = 8  # Distance (in pixel) from which the selection of curve occurs


class AutomaticAnalysis(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        """
        Creates the specific windows in the main parent window and takes
        the figure, canvas, axe and toolbar to be used.
        """
        tk.Frame.__init__(self, parent, *args, **kwargs)  # Frame initialization
        self.parent = parent  # Saving the parent's pointer
        self.parent.wm_title('Spectral Analysis | AUTOMATIC ANALYSIS')  # Changing the parent's title

        self.click = False  # Check variable if button is pressed
        self.grange = {}  # Graphical (plot) range (by area name)
        self.grangeid = None
        self.poptsave = {}  # Peaks parameters (by name)
        self.editoff = False  # Check variable if edition of parameters is allowed
        self.tabId = None  # Name of the focused tab
        self.peaks = []  # Axes of the shown peaks

        # Id of the closest value in the x datas from the given area boundaries
        self. areaBoundsId = {}
        for area in self.parent.areaBounds:
            x1id = abs(self.parent.data[0] - self.parent.areaBounds[area][0]).argmin()
            x2id = abs(self.parent.data[0] - self.parent.areaBounds[area][1]).argmin()
            self.areaBoundsId[area] = sorted([x1id, x2id])

        # Interface
        self.makeInterface()

    def makeInterface(self):
        # Deconvolute button
        self.refreshbutton = tk.Button(
            master=self, text='Deconvolute',
            command=self.refreshCallback)
        self.refreshbutton.grid(row=0, column=0, pady=10, padx=5)

        self.parent.bind("<Return>", lambda event: self.apply())

        # Save Button
        self.savebutton = tk.Button(
            master=self, text='Save Param.',
            command=self.saveCallback)
        self.savebutton.grid(row=2, column=0, pady=10, padx=5)

        # Event connection
        self.eventId = []
        self.eventId.append(self.parent.canvas.mpl_connect('button_press_event', self.on_click))
        self.eventId.append(self.parent.canvas.mpl_connect('button_release_event', self.on_release))
        self.eventId.append(self.parent.canvas.mpl_connect('pick_event', self.on_pick))
        self.eventId.append(self.parent.canvas.mpl_connect('motion_notify_event', self.on_motion))

        # Colors for area boundaries
        color = ['r', 'g', 'm', 'c', 'y']

        i = -1  # Iterator for color
        for area in self.parent.areaBounds:
            i += 1
            self.grange[area] = [
                self.parent.ax.plot(
                    [min(self.parent.areaBounds[area]), min(self.parent.areaBounds[area])],
                    [min(self.parent.data[1]), max(self.parent.data[1])],
                    "-" + color[i%len(color)], picker=PICKERSIZE, label=area)[0],
                self.parent.ax.plot(
                    [max(self.parent.areaBounds[area]), max(self.parent.areaBounds[area])],
                    [min(self.parent.data[1]), max(self.parent.data[1])],
                    "-" + color[i%len(color)], picker=PICKERSIZE)[0],
                self.parent.ax.plot(
                    [min(self.parent.areaBounds[area]), max(self.parent.areaBounds[area])],
                    [self.parent.data[1][abs(self.parent.data[0] - min(self.parent.areaBounds[area])).argmin()],
                     self.parent.data[1][abs(self.parent.data[0] - max(self.parent.areaBounds[area])).argmin()]],
                    "-" + color[i%len(color)])[0]]

        plt.legend()

        # Notebook creation
        self.notebook = Notebook(self)
        self.notebook.grid(column=1, row=0, rowspan=100, sticky=tk.NW+tk.SE)

        # Notebook's tabs creation
        self.tabs = {}
        for area in self.parent.areaBounds:
            self.tabs[area] = tk.Frame(master=self.notebook)
            self.notebook.add(self.tabs[area], text=area)

        # Defining the tab in focus
        self.parammod()

        # Linking the change in tab focus with the callback
        self.notebook.bind('<<NotebookTabChanged>>', self.parammod)

        # Add/remove peak button creation
        self.addPeak = tk.Button(master=self, text='Add Peak', command=self.add)
        self.addPeak.grid(column=2, row=0, rowspan=100, sticky=tk.NW+tk.SE)

        self.addPeak = tk.Button(master=self, text='Remove Peak', command=self.remove)
        self.addPeak.grid(column=3, row=0, rowspan=100, sticky=tk.NW+tk.SE)

    def add(self):
        """ Add peak to the focused tab """
        area = self.tabId
        rndWaveNb = uniform(self.parent.areaBounds[area][0], self.parent.areaBounds[area][1])
        self.parent.guess[area] = np.concatenate((self.parent.guess[area], np.array([0.0001, 5, rndWaveNb])))
        self.parent.peakBounds[area][0] = np.concatenate((self.parent.peakBounds[area][0], np.array([0.0001, 0, -np.inf])))
        self.parent.peakBounds[area][1] = np.concatenate((self.parent.peakBounds[area][1], np.array([np.inf, np.inf, np.inf])))
        self.parammod()

    def remove(self):
        """ Remove last peak from the focused tab """
        area = self.tabId
        if len(self.parent.guess[area]) > 3:
            self.parent.guess[area] = self.parent.guess[area][:-3]
            self.parent.peakBounds[area][0] = self.parent.peakBounds[area][0][:-3]
            self.parent.peakBounds[area][1] = self.parent.peakBounds[area][1][:-3]
            self.parammod()

    def refreshCallback(self):
        """ Automatic deconvolution callback """
        for i in range(len(self.peaks)):
            self.peaks[i].remove()
        self.peaks = []

        guess = deepcopy(self.parent.guess)
        peakBounds = deepcopy(self.parent.peakBounds)

        for area in self.parent.areaBounds:  # For every area
            # Restraining the ROI to the analysed area
            # Getting Id of the boundaries
            x = self.parent.data[0][
                self.areaBoundsId[area][0]:
                self.areaBoundsId[area][1]+1]
            y = self.parent.data[1][
                self.areaBoundsId[area][0]:
                self.areaBoundsId[area][1]+1]

            # Calculating and subtracting the baseline of the ROI
            bl = sf.baseline(x, y, method=self.parent.blMethod)
            y = y - bl

            # Deconvolution

            # Bone specific
            if area == 'amide':
                # Add residual MMA peak if still present
                arg1730x = abs(x - 1730).argmin()
                arg1760x = abs(x - 1760).argmin()
                yp = sf.derivate(x[arg1730x: arg1760x + 1], y[arg1730x: arg1760x + 1])
                yp_pos = deepcopy(yp)
                yp_pos[yp < 0] = 0
                yp_neg = deepcopy(yp)
                yp_neg[yp > 0] = 0

                if np.trapz(yp_pos, x[arg1730x: arg1760x + 1]) > 0.5 * abs(np.trapz(yp_neg, x[arg1730x: arg1760x + 1])):
                    print("Adding MMA residual peak")
                    guess['amide'] = np.append(guess['amide'], [0.03, 5, 1745])
                    peakBounds['amide'][0] = np.append(peakBounds['amide'][0], [0.0001, 0.0001, 1735])
                    peakBounds['amide'][1] = np.append(peakBounds['amide'][1], [np.inf, 15, 1755])

            # First guess of the height of peaks at 90% intensity of curve
            for wavenb in [guess[area][i] for i in range(len(guess[area])) if (i-2) % 3 == 0]:
                # Height constraint of the near 1633 peak
                try:
                    wavenbs = np.array([guess[area][i * 3 + 2] for i in
                                        range(len(guess[area]) // 3)])
                    arg1633 = abs(wavenbs - wavenb).argmin()
                    arg1633x = abs(x - wavenb).argmin()
                    guess[area][arg1633 * 3] = y[arg1633x] * 0.9
                except Exception as e:
                    print(str(e))

                """
                # Height constraint of the near 1600 peak
                try:
                    heights = np.array([self.parent.peakBounds['amide'][1][i*3] for i in range(int(len(self.parent.peakBounds['amide'][1])/3))])
                    wavenbs = np.array([self.parent.peakBounds['amide'][1][i*3+2] for i in range(int(len(self.parent.peakBounds['amide'][1])/3))])
                    arg1600 = abs(wavenbs - 1600).argmin()
                    arg1600x = abs(x - 1600).argmin()
                    self.parent.peakBounds['amide'][1][arg1600*3] = y[arg1600x]*2./3
                    self.parent.guess['amide'][arg1600*3] = y[arg1600x]*2./3*0.9
                except Exception as e:
                    print(str(e))
                    
                # Height constraint of the near 1633 peak
                try:
                    heights = np.array([self.parent.peakBounds['amide'][1][i * 3] for i in
                                        range(len(self.parent.peakBounds['amide'][1]) // 3)])
                    wavenbs = np.array([self.parent.peakBounds['amide'][1][i * 3 + 2] for i in
                                        range(len(self.parent.peakBounds['amide'][1]) // 3)])
                    arg1633 = abs(wavenbs - 1633).argmin()
                    arg1633x = abs(x - 1633).argmin()
                    self.parent.peakBounds['amide'][0][arg1633 * 3] = y[arg1633x] * 0.5
                except Exception as e:
                    print(str(e))
                """
            # if area == 'v1v3PO4':
            #     print(peakBounds[area][1][-3])
            #     # Height constraint of the near 1060 peak (based on the 1030 peak)
            #     try:
            #         heights = np.array([peakBounds['v1v3PO4'][1][i*3] for i in range(int(len(peakBounds['v1v3PO4'][1])/3))])
            #         wavenbs = np.array([peakBounds['v1v3PO4'][1][i*3+2] for i in range(int(len(peakBounds['v1v3PO4'][1])/3))])
            #         arg1030x = abs(x - 1030).argmin()
            #         arg1060 = abs(wavenbs - 1060).argmin()
            #         peakBounds['v1v3PO4'][1][arg1060*3] = y[arg1030x]*1./3
            #
            #         arg1030 = abs(wavenbs - 1030).argmin()
            #         arg1110 = abs(wavenbs - 1110).argmin()
            #         arg1110x = abs(x - 1110).argmin()
            #         peakBounds['v1v3PO4'][0][arg1030*3] = y[arg1030x]*4./5
            #         peakBounds['v1v3PO4'][0][arg1110*3] = y[arg1110x]*4./5
            #         guess['v1v3PO4'][arg1030*3] = y[arg1030x]*4./5
            #         guess['v1v3PO4'][arg1110*3] = y[arg1110x]*4./5
            #
            #         # heights = np.array([self.parent.peakBounds['v1v3PO4'][1][i*3] for i in range(int(len(self.parent.peakBounds['v1v3PO4'][1])/3))])
            #         # wavenbs = np.array([self.parent.peakBounds['v1v3PO4'][1][i*3+2] for i in range(int(len(self.parent.peakBounds['v1v3PO4'][1])/3))])
            #         # arg1030x = abs(x - 1030).argmin()
            #         # arg1060 = abs(wavenbs - 1060).argmin()
            #         #
            #         # arg1030 = abs(wavenbs - 1030).argmin()
            #         # arg1110 = abs(wavenbs - 1110).argmin()
            #         # arg1110x = abs(x - 1110).argmin()
            #         # self.parent.peakBounds['v1v3PO4'][1][arg1060*3] = y[arg1110x]*3./4
            #         # self.parent.peakBounds['v1v3PO4'][0][arg1030*3] = y[arg1030x]*4./5
            #         # self.parent.peakBounds['v1v3PO4'][0][arg1110*3] = y[arg1110x]*4./5
            #         # self.parent.guess['v1v3PO4'][arg1030*3] = y[arg1030x]*4./5
            #         # self.parent.guess['v1v3PO4'][arg1110*3] = y[arg1110x]*4./5
            #
            #     except Exception as e:
            #         print(str(e))

            popt = sf.spectrumFit(self.parent.function,
                                  x, y, guess=guess[area],
                                  bounds=peakBounds[area])

            # Plotting the peaks
            color = ['#1dba27', '#dc9100', '#9d27d2', '#1a5bb4', '#ddd526']
            marker = ['-', '--']
            for j in range(0, int(len(popt) / 3 + 1)):
                self.peaks.append(self.parent.ax.plot(x, self.parent.function(x, *popt[j * 3 - 3:j * 3]) + bl, marker[(j-1)//5 % 2], color=color[j % 5])[0])
                self.peaks.append(self.parent.ax.plot(x, self.parent.function(x, *popt) + bl, 'k:', linewidth=3)[0])

            # Saving the ROI peaks
            self.poptsave[area] = popt

        if self.parent.display:
            self.parent.canvas.draw()  # Updating the canvas

    def on_click(self, event):
        """ Mouse button clicked callback """
        if self.parent.toolbar._active is None:  # If no tools from toolbar selected
            self.click = True  # Activate the pressed state

    def on_motion(self, event):
        """
        Mouse on motion callback.
        Drags and updates the selected area boundary.
        """
        # If pressed state is active and no tools from toolbar is selected
        if self.click is True and self.parent.toolbar._active is None:
            x = event.xdata  # Getting the x coordinate of the mouse
            # Finding the nearest value in the spectrum's x axis
            xid = abs(self.parent.data[0] - x).argmin()
            # Updates the visualization of the area boundaries
            self.grange[self.grangeid][self.grangelim].set_data(
                [self.parent.data[0][xid], self.parent.data[0][xid]],
                [min(self.parent.data[1]), max(self.parent.data[1])])
            # Updates the value of the boundary
            if self.grangelim == 0:
                self.parent.areaBounds[self.grangeid][0] = x
                self.areaBoundsId[self.grangeid][0] = xid
            else:
                self.parent.areaBounds[self.grangeid][1] = x
                self.areaBoundsId[self.grangeid][1] = xid
            # Updates the linking line of the area boundaries
            self.grange[self.grangeid][2].set_data(
                [self.parent.areaBounds[self.grangeid][0], self.parent.areaBounds[self.grangeid][1]],
                self.parent.data[1][[self.areaBoundsId[self.grangeid][0],
                                       self.areaBoundsId[self.grangeid][1]]])

            if self.parent.display:
                self.parent.canvas.draw()  # Canvas update

    def on_release(self, event):
        """ Deactivate the pressed state when mouse button is released """
        self.click = False

    def on_pick(self, event):
        """ Line selection callback """
        # If no tools from toolbar is selected
        if self.parent.toolbar._active is None:
            # Thickness of previously selected boundary to default
            if self.grangeid is not None:
                for line in self.grange[self.grangeid]:
                    plt.setp(line, linewidth=1)

            # Finding the selected line index
            thisline = event.artist
            for area in self.grange:
                if thisline is self.grange[area][0]:
                    self.grangeid = area
                    self.grangelim = 0
                    break
                if thisline is self.grange[area][1]:
                    self.grangeid = area
                    self.grangelim = 1
                    break

            # Thicker selected line
            for line in self.grange[self.grangeid]:
                plt.setp(line, linewidth=2)

            if self.parent.display:
                self.parent.canvas.draw()  # Canvas update

    def parammod(self, *args):
        """ Defines the tab in focus """
        # Finding the tab in focus
        for area in self.parent.areaBounds:
            if self.notebook.select() == str(self.tabs[area]):
                self.tabId = area
                master = self.tabs[area]
                break
        # Grid window
        self.subwin = tk.Frame(master=master)
        self.subwin.grid(
            sticky=tk.W + tk.S + tk.N + tk.E, row=0, column=0, padx=0)

        # Labels
        self.heightlab = tk.Label(master=self.subwin, text='Max Height')
        self.heightlab.grid(column=0, row=1, sticky=tk.E, padx=2, pady=2)

        self.fwhmlab = tk.Label(master=self.subwin, text='Max FWHM/2')
        self.fwhmlab.grid(column=0, row=2, sticky=tk.E, padx=2, pady=2)

        self.wavenblab = tk.Label(master=self.subwin, text='Wave Number')
        self.wavenblab.grid(column=0, row=3, sticky=tk.E, padx=2, pady=2)

        self.wavenblab = tk.Label(
            master=self.subwin, text='Wave Nb range (+/-)')
        self.wavenblab.grid(column=0, row=4, sticky=tk.E, padx=2, pady=2)

        # Entries (Defining a visual grid of parameters: each column = a peak)
        # tabId defines the current focused tab
        self.guess = self.parent.guess[self.tabId]
        self.bounds = self.parent.peakBounds[self.tabId]
        nbPeak = int(len(self.guess) / 3)

        self.heightsel = [tk.StringVar(value='', ) for i in range(nbPeak)]
        self.fwhmsel = [tk.StringVar(value='') for i in range(nbPeak)]
        self.wavenbsel = [tk.StringVar(value='') for i in range(nbPeak)]
        self.wavenbrangesel = [tk.StringVar(value='') for i in range(nbPeak)]

        self.heightedit = [tk.Entry(master=self.subwin,
                                    width=4,
                                    textvariable=self.heightsel[i])
                           for i in range(nbPeak)]
        [self.heightedit[i].grid(column=i + 1, row=1, padx=2, pady=2)
            for i in range(nbPeak)]

        self.fwhmedit = [tk.Entry(master=self.subwin,
                                  width=4,
                                  textvariable=self.fwhmsel[i])
                         for i in range(nbPeak)]
        [self.fwhmedit[i].grid(column=i + 1, row=2, padx=2, pady=2)
            for i in range(nbPeak)]

        self.wavenbedit = [tk.Entry(master=self.subwin,
                                    width=4,
                                    textvariable=self.wavenbsel[i])
                           for i in range(nbPeak)]
        [self.wavenbedit[i].grid(column=i + 1, row=3, padx=2, pady=2)
            for i in range(nbPeak)]

        self.wavenbrangeedit = [tk.Entry(master=self.subwin,
                                         width=4,
                                         textvariable=self.wavenbrangesel[i])
                                for i in range(nbPeak)]
        [self.wavenbrangeedit[i].grid(column=i + 1, row=4, padx=2, pady=2)
            for i in range(nbPeak)]

        self.editoff = True
        [self.heightsel[i].set(str(self.bounds[1][i * 3]))
            for i in range(nbPeak)]
        [self.fwhmsel[i].set(str(self.bounds[1][i * 3 + 1]))
            for i in range(nbPeak)]
        [self.wavenbsel[i].set(str(self.guess[i * 3 + 2]))
            for i in range(nbPeak)]
        [self.wavenbrangesel[i].set(
            str((self.bounds[1][i * 3 + 2] - self.bounds[0][i * 3 + 2]) / 2))
            for i in range(nbPeak)]
        self.editoff = False

        # Connecting the entries events (partial is used to send i, the peak
        # index, to the event callback)
        [self.heightsel[i].trace("w", partial(self.heightselCallback, i))
            for i in range(nbPeak)]
        [self.fwhmsel[i].trace("w", partial(self.fwhmselCallback, i))
            for i in range(nbPeak)]
        [self.wavenbsel[i].trace("w", partial(self.wavenbselCallback, i))
            for i in range(nbPeak)]
        [self.wavenbrangesel[i].trace(
            "w", partial(self.wavenbrangeselCallback, i))
            for i in range(nbPeak)]

    # Entries callbacks
    def heightselCallback(self, i, *args):
        """
        Modify the concerned parameter of the peak i (from the focused tab)
        with the content of the entry
        """
        self.parent.peakBounds[self.tabId][1][i * 3] = float(self.heightsel[i].get())
        return

    def fwhmselCallback(self, i, *args):
        """
        Modify the concerned parameter of the peak i (from the focused tab)
        with the content of the entry
        """
        self.parent.peakBounds[self.tabId][1][i * 3 + 1] = float(self.fwhmsel[i].get())
        return

    def wavenbselCallback(self, i, *args):
        """
        Modify the concerned parameter of the peak i (from the focused tab)
        with the content of the entry
        """
        self.parent.guess[self.tabId][i * 3 + 2] = float(self.wavenbsel[i].get())
        self.parent.peakBounds[self.tabId][0][i * 3 + 2] = \
            self.parent.guess[self.tabId][i * 3 + 2] - \
            float(self.wavenbrangesel[i].get())
        self.parent.peakBounds[self.tabId][1][i * 3 + 2] = \
            self.parent.guess[self.tabId][i * 3 + 2] + \
            float(self.wavenbrangesel[i].get())
        return

    def wavenbrangeselCallback(self, i, *args):
        """
        Modify the concerned parameter of the peak i (from the focused tab)
        with the content of the entry
        """
        self.parent.peakBounds[self.tabId][0][i * 3 + 2] = \
            self.parent.guess[self.tabId][i * 3 + 2] - \
            float(self.wavenbrangesel[i].get())
        self.parent.peakBounds[self.tabId][1][i * 3 + 2] = \
            self.parent.guess[self.tabId][i * 3 + 2] + \
            float(self.wavenbrangesel[i].get())
        return

    def saveCallback(self, *args):
        """
        Compute and save the parameters used for the FTIRM analysis of
        bone samples ('Mineralization Index', 'Carbonation', 'Crystallinity',
        'Mineral Maturity', 'Crystal Maturity').
        The output filename is "'originalDataFileName'_results.csv" or
        "'originalDataFileName'_results(i).csv" with i a number starting from 1
        if the file already exists.
        """
        # v1v3PO4 parameters
        for area in self.parent.areaBounds:
            if area == 'v1v3PO4':
                # Restraining datas to the ROI
                x = self.parent.data[0][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                y = self.parent.data[1][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                bl = sf.baseline(x, y, method=self.parent.blMethod)
                y = y - bl

                # Calculating the area/integral of the concerned peaks
                # Resampling the datas for better area calculation

                v1v3PO4_area = abs(np.trapz(y, x))
                if self.parent.limReached is True:
                    v1v3PO4_area *= (1 - np.exp(0.03878*self.parent.lenLim - 6.652))

                wavenbs = np.array([self.poptsave[area][j * 3 + 2] for j in range(int(len(self.poptsave[area]) / 3))])

                peak1030id_intensity = abs(x - 1030).argmin()
                peak1110id_intensity = abs(x - 1110).argmin()

                peak1030_intensity = y[peak1030id_intensity]
                peak1110_intensity = y[peak1110id_intensity]

                peak1030id = abs(wavenbs - 1030).argmin()
                peak1110id = abs(wavenbs - 1110).argmin()
                print('1110 cm-1 peak found at : %f cm-1'%(self.poptsave[area][peak1110id * 3 + 2]))
                print('1030 cm-1 peak found at : %f cm-1'%(self.poptsave[area][peak1030id * 3 + 2]))
                x = np.linspace(-500, 500, 10000)
                peak1030_area = abs(
                    np.trapz(
                        self.parent.function(
                            x, self.poptsave[area][peak1030id * 3], self.poptsave[area][peak1030id * 3 + 1], 0),
                        x))
                peak1110_area = abs(
                    np.trapz(
                        self.parent.function(
                            x, self.poptsave[area][peak1110id * 3], self.poptsave[area][peak1110id * 3 + 1], 0),
                        x))


                #print peak1030_area
                #print self.poptsave[i][peak1030id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1030id * 3+1]/np.sqrt(np.log(2)))**2)
                #print peak1110_area
                #print self.poptsave[i][peak1110id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1110id * 3+1]/np.sqrt(np.log(2)))**2)

        # v4PO4 parameters
            if area == 'v4PO4':
                # Restraining datas to the ROI
                x = self.parent.data[0][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                y = self.parent.data[1][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                bl = sf.baseline(x, y, method=self.parent.blMethod)
                y = y - bl

                # Calculating the fwhm of the concerned peaks
                wavenbs = np.array([self.poptsave[area][j * 3 + 2] for j in range(int(len(self.poptsave[area]) / 3))])

                peak604id = abs(wavenbs - 604).argmin()
                xtest = np.linspace(0, 100, 2000)
                peak604 = np.array(sf.gaussianSum(xtest, self.poptsave[area][peak604id * 3],
                                                  self.poptsave[area][peak604id * 3 + 1], 0))
                peak604_fwhm = xtest[abs(peak604 - max(peak604)/2).argmin()] * 2

        # CO3 parameters
            if area == 'CO3':
                # Restraining datas to the ROI
                x = self.parent.data[0][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                y = self.parent.data[1][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                yorig = y
                xorig = x

                idmin = len(y) - 1
                y = yorig[:idmin+1]
                x = xorig[:idmin+1]
                idstart = 0
                a = (yorig[idmin] - yorig[idstart])/(xorig[idmin] - xorig[idstart])
                b = yorig[0] - a*xorig[0]
                y = y - sf.baseline(x, y)

                """xstand = x
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
                        break"""

                # Calculating the area/integral of the concerned peaks
                CO3_area = np.trapz(y, x)

        # amide parameters
            if area == 'amide':
                # Restraining datas to the ROI
                x = self.parent.data[0][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                y = self.parent.data[1][self.areaBoundsId[area][0]:
                                        self.areaBoundsId[area][1]+1]
                bl = sf.baseline(x, y, method=self.parent.blMethod)
                y = y - bl

                # Calculating the area/integral of the concerned peaks
                wavenbs = np.array([self.poptsave[area][j * 3 + 2] for j in range(int(len(self.poptsave[area]) / 3))])

                peak1660id = abs(wavenbs - 1660).argmin()
                peak1690id = abs(wavenbs - 1690).argmin()

                peak1660id_intensity = abs(x - 1660).argmin()
                peak1690id_intensity = abs(x - 1690).argmin()

                peak1660_intensity = y[peak1660id_intensity]
                peak1690_intensity = y[peak1690id_intensity]

                x = np.linspace(-500, 500, 10000)
                peak1660_area = abs(
                    np.trapz(
                        self.parent.function(
                            x, self.poptsave[area][peak1660id * 3], self.poptsave[area][peak1660id * 3 + 1], 0),
                        x))
                peak1690_area = abs(
                    np.trapz(
                        self.parent.function(
                            x, self.poptsave[area][peak1690id * 3], self.poptsave[area][peak1690id * 3 + 1], 0),
                        x))
                #print peak1690_area
                #print self.poptsave[i][peak1690id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1690id * 3+1]/np.sqrt(np.log(2)))**2)
                #print peak1660_area
                #print self.poptsave[i][peak1660id * 3]*np.sqrt(np.pi*(self.poptsave[i][peak1660id * 3+1]/np.sqrt(np.log(2)))**2)

        # Amide I area
        try:
            x1 = 1592
            x2 = 1730
            x1id = abs(self.parent.data[0] - x1).argmin()
            x2id = abs(self.parent.data[0] - x2).argmin()
            x = self.parent.data[0][min(x1id, x2id):max(x1id, x2id)+1]
            y = self.parent.data[1][min(x1id, x2id):max(x1id, x2id)+1]
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
            MI_intensity = peak1030_intensity / peak1660_intensity
            if MI < 1.96 or MI > 8.56:
                failure = True
        except:
            MI = None
            MI_intensity = None
            pass
        try:
            carbonation = CO3_area / v1v3PO4_area
            if carbonation < 0.00301 or carbonation > 0.0120:
                failure = True
        except:
            carbonation = None
            pass
        try:
            crystallinity = 1./peak604_fwhm
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
            CM_intensity = peak1660_intensity / peak1690_intensity
            if CM < 1.83 or CM > 6.58:
                failure = True
        except:
            CM = None
            CM_intensity = None
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
                 'Crystallinity', 'Mineral Maturity', 'Collagen Maturity', 'Collagen Maturity (intensity)', 'Mineralization Index (intensity)', 'Mineral Maturity (intensity)', 'v1v3PO4 area', 'Cut length']
        # Second line
        results = [MI, carbonation, crystallinity, MM, CM, CM_intensity, MI_intensity, MM_intensity, v1v3PO4_area, self.parent.lenLim]
        # Extracting the original filname
        path = os.path.dirname(self.parent.filepath)
        name = os.path.splitext(os.path.basename(self.parent.filepath))[0]
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
        xmin = min([self.parent.areaBounds[area][0] for area in self.parent.areaBounds])
        xmax = max([self.parent.areaBounds[area][1] for area in self.parent.areaBounds])
        totRange = xmax-xmin
        xmin -= 0.02*totRange
        xmax += 0.02*totRange

        origxlim = self.parent.ax.get_xlim()
        origsize = plt.gcf().get_size_inches()
        self.parent.ax.set_xlim(xmax, xmin)
        plt.gcf().set_size_inches(origsize[0]*6, origsize[1]*2)
        #if not failure:
        plt.savefig(path + '/Results/' + name + '_results_spectrum.png', dpi=300)
        #else:
            #plt.savefig(path + '/Results_fail/' + name + '_results_spectrum.png', dpi=300)
        plt.gcf().set_size_inches(origsize[0], origsize[1])
        self.parent.ax.set_xlim(origxlim)

    def _quit(self):
        """ Properly quits the automatic analysis tool """
        self.parent.wm_title('Spectral Analysis')
        for i in self.eventId:
            self.parent.canvas.mpl_disconnect(i)
        self.parent.unbind("<Return>")
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.destroy()
        return
