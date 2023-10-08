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
Manual Analysis class.
Defines a manual analysis window and tools.
"""

import numpy as np
import matplotlib.pyplot as plt
import spectrumFit as sf
import tkinter as tk


FILENAME = './configure.ini'
pickersize = 8  # Distance (in pixel) from which the selection of curve occurs


class ManualAnalysis(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        """
        Creates the specific windows in the main parent window and takes
        the figure, canvas, axe and toolbar to be used.
        """
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.parent.wm_title('Spectral Analysis | MANUAL ANALYSIS')

        self.dataROI = [self.parent.data[0][self.parent.bounds[0]: self.parent.bounds[1]],
                        self.parent.data[1][self.parent.bounds[0]: self.parent.bounds[1]]]

        self.nbpeak = 0  # Total number of peaks
        self.lineid = None

        self.height = np.array([])  # Stores the height of each peak
        self.fwhm = np.array([])  # Stores the fwhm of each peak
        self.wavenb = np.array([])  # Stores the wave number of each peak
        self.heightrange = np.array([])  # Stores the height range of each peak
        self.fwhmrange = np.array([])  # Stores the fwhm range of each peak
        self.wavenbrange = np.array([])  # Stores wave nb range of each peak

        self.function = self.parent.function

        # Text value of the selected peak's parameters
        self.heightsel = tk.StringVar(value='')
        self.fwhmsel = tk.StringVar(value='')
        self.wavenbsel = tk.StringVar(value='')
        self.heightrangesel = tk.StringVar(value='')
        self.fwhmrangesel = tk.StringVar(value='')
        self.wavenbrangesel = tk.StringVar(value='')

        self.refreshcount = 0  # Check variable if deconvolution button pressed
        self.peakcreation = False  # Check variable if peak is creating

        # Grid window (nbpeak/selected peak/refresh)
        self.subwin = [
            tk.Frame(master=self) for i in [tk.LEFT, None, tk.RIGHT]]
        [self.subwin[i].grid(
            sticky=tk.W + tk.S + tk.N + tk.E, row=0, column=i, padx=20)
            for i in range(0, len(self.subwin))]
        [self.columnconfigure(i, weight=1)
            for i in range(0, len(self.subwin))]

        # Labels
        self.peaklab = tk.Label(
            master=self.subwin[0],
            text='Number of Peaks: {}'.format(self.nbpeak))
        self.peaklab.pack(fill=tk.BOTH, expand=1)

        self.selectedPeaklab = tk.Label(
            master=self.subwin[1], text='Selected Peak')
        self.selectedPeaklab.grid(
            column=1, row=0, columnspan=1, sticky=tk.N, pady=10)

        self.rangelab = tk.Label(master=self.subwin[1], text='Range')
        self.rangelab.grid(column=3, row=0, columnspan=1, sticky=tk.N, pady=10)

        self.heightlab = tk.Label(master=self.subwin[1], text='Height')
        self.heightlab.grid(column=0, row=1, sticky=tk.E, padx=3, pady=3)

        self.fwhmlab = tk.Label(master=self.subwin[1], text='FWHM')
        self.fwhmlab.grid(column=0, row=2, sticky=tk.E, padx=3, pady=3)

        self.wavenblab = tk.Label(master=self.subwin[1], text='Wave Number')
        self.wavenblab.grid(column=0, row=3, sticky=tk.E, padx=3, pady=3)

        self.pom0lab = tk.Label(master=self.subwin[1], text='+/-')
        self.pom0lab.grid(column=2, row=1, padx=3, pady=3)

        self.pom1lab = tk.Label(master=self.subwin[1], text='+/-')
        self.pom1lab.grid(column=2, row=2, padx=3, pady=3)

        self.pom2lab = tk.Label(master=self.subwin[1], text='+/-')
        self.pom2lab.grid(column=2, row=3, padx=3, pady=3)

        # Entries
        self.heightedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.heightsel)
        self.heightedit.grid(column=1, row=1, padx=3, pady=3)

        self.fwhmedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.fwhmsel)
        self.fwhmedit.grid(column=1, row=2, padx=3, pady=3)

        self.wavenbedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.wavenbsel)
        self.wavenbedit.grid(column=1, row=3, padx=3, pady=3)

        self.heightrangeedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.heightrangesel)
        self.heightrangeedit.grid(column=3, row=1, padx=3, pady=3)

        self.fwhmrangeedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.fwhmrangesel)
        self.fwhmrangeedit.grid(column=3, row=2, padx=3, pady=3)

        self.wavenbrangeedit = tk.Entry(
            master=self.subwin[1], width=5, textvariable=self.wavenbrangesel)
        self.wavenbrangeedit.grid(column=3, row=3, padx=3, pady=3)

        # Entries/StringVar event connection
        self.heightsel.trace("w", self.heightselCallback)
        self.fwhmsel.trace("w", self.fwhmselCallback)
        self.wavenbsel.trace("w", self.wavenbselCallback)
        self.heightrangesel.trace("w", self.heightrangeselCallback)
        self.fwhmrangesel.trace("w", self.fwhmrangeselCallback)
        self.wavenbrangesel.trace("w", self.wavenbrangeselCallback)

        # Buttons
        self.deletebutton = tk.Button(
            master=self.subwin[2], text='Delete selected peak',
            command=self.deleteCallback)
        self.deletebutton.pack(side=tk.TOP, anchor=tk.W)

        self.refreshbutton = tk.Button(
            master=self.subwin[2], text='Deconvolve',
            command=self.refreshCallback)
        self.refreshbutton.pack(side=tk.TOP, anchor=tk.W, pady=20)

        self.parent.bind("<Return>", lambda event: self.apply())

        # Canvas events connactions
        self.eventId = []
        self.eventId.append(self.parent.canvas.mpl_connect('button_press_event', self.on_click))
        self.eventId.append(self.parent.canvas.mpl_connect('pick_event', self.on_pick))

    def heightselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.height[self.lineid] = float(self.heightsel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def fwhmselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.fwhm[self.lineid] = float(self.fwhmsel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def wavenbselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.wavenb[self.lineid] = float(self.wavenbsel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def heightrangeselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.heightrange[self.lineid] = float(self.heightrangesel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def fwhmrangeselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.fwhmrange[self.lineid] = float(self.fwhmrangesel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def wavenbrangeselCallback(self, *args):
        """ Change value by the one in the entry widget and update datas """
        self.wavenbrange[self.lineid] = float(self.wavenbrangesel.get())
        plt.gca().get_lines()[self.lineid + 1].set_ydata(
            self.function(self.dataROI[0], self.height[self.lineid],
                           self.fwhm[self.lineid] / 2,
                           self.wavenb[self.lineid]))
        if self.peakcreation is False:
            self.parent.canvas.draw()

    def deleteCallback(self):
        """ Deletes the selected peak """
        if self.lineid is None:  # If no peak selected
            return
        self.nbpeak -= 1  # Decrease total number of peaks
        self.peaklab.config(text='Number of Peaks : {}'.format(self.nbpeak))
        # Updates the value

        plt.gca().get_lines().pop(self.lineid + 1).remove()
        # Removes the plot of the selected peak

        # Deletes the parameters's value of the selected peak
        self.height = np.delete(self.height, self.lineid, None)
        self.wavenb = np.delete(self.wavenb, self.lineid, None)
        self.fwhm = np.delete(self.fwhm, self.lineid, None)
        self.heightrange = np.delete(self.heightrange, self.lineid, None)
        self.wavenbrange = np.delete(self.wavenbrange, self.lineid, None)
        self.fwhmrange = np.delete(self.fwhmrange, self.lineid, None)

        self.peakcreation = True  # Activation of the peake creation state
        # (deactivates the automatic update of the values due to the entry
        # text change)
        self.heightsel.set('')
        self.wavenbsel.set('')
        self.fwhmsel.set('')
        self.heightrangesel.set('')
        self.wavenbrangesel.set('')
        self.fwhmrangesel.set('')
        self.peakcreation = False  # Deactivation of peak creation state

        self.lineid = None  # No line selected
        self.parent.canvas.draw()  # Update canvas

    def refreshCallback(self):
        """ Deconvolve the spectrum with the created peaks """
        guess = np.zeros(3 * len(self.height))  # Initial guess of the params
        bound_inf = np.zeros(3 * len(self.height))  # Inf boundary
        bound_sup = np.zeros(3 * len(self.height))  # Sup boundary

        # Getting the values for each peak
        for i in range(len(self.height)):
            guess[i * 3] = self.height[i]
            guess[i * 3 + 1] = self.fwhm[i] / 2
            guess[i * 3 + 2] = self.wavenb[i]
            bound_inf[i * 3] = self.height[i] - self.heightrange[i]
            bound_inf[i * 3 + 1] = (self.fwhm[i] - self.fwhmrange[i]) / 2
            bound_inf[i * 3 + 2] = self.wavenb[i] - self.wavenbrange[i]
            bound_sup[i * 3] = self.height[i] + self.heightrange[i]
            bound_sup[i * 3 + 1] = (self.fwhm[i] + self.fwhmrange[i]) / 2
            bound_sup[i * 3 + 2] = self.wavenb[i] + self.wavenbrange[i]
        bounds = [bound_inf, bound_sup]

        # Curve fitting
        popt = sf.spectrumFit(
            self.function, self.dataROI[0], self.dataROI[1],
            guess=guess, bounds=bounds)

        # Updates the results for each peak
        for i in range(len(self.height)):
            self.height[i] = popt[i * 3]
            self.fwhm[i] = popt[i * 3 + 1] * 2
            self.wavenb[i] = popt[i * 3 + 2]

        # Updates values for selected peak
        if hasattr(self, 'lineid'):
            if self.lineid is not None:
                self.heightsel.set(str(self.height[self.lineid]))
                self.wavenbsel.set(str(self.wavenb[self.lineid]))
                self.fwhmsel.set(str(self.fwhm[self.lineid]))
                self.heightrangesel.set(str(self.heightrange[self.lineid]))
                self.wavenbrangesel.set(str(self.wavenbrange[self.lineid]))
                self.fwhmrangesel.set(str(self.fwhmrange[self.lineid]))

        # Updates each peak
        lines = plt.gca().get_lines()
        for i in range(1, len(lines)):
            lines[i].set_ydata(
                self.function(self.dataROI[0], *popt[i * 3 - 3:i * 3]))
        if self.refreshcount != 0:  # If deconv. done before, update the fit
            plt.gca().get_lines()[-1].set_ydata(
                self.function(self.dataROI[0], *popt))
        else:  # If first deconvolution, creates the fit (sum of each peak)
            self.parent.ax.plot(
                self.dataROI[0],
                self.function(self.dataROI[0], *popt),
                linewidth=1.0, picker=pickersize)

        self.refreshcount += 1  # Increments the refresh count
        self.parent.canvas.draw()  # Updates the canvas

    # Canvas event
    def on_click(self, event):
        """
        Mouse click callback.
        Right click = creation.
        """
        # Creation of a peak
        # Right button, no tools active
        if event.inaxes is not None and event.button == 3 \
                and self.parent.toolbar._active is None:
            if event.xdata > max(self.dataROI[0]):
                event.xdata = max(self.dataROI[0])
            if event.xdata < min(self.dataROI[0]):
                event.xdata = min(self.dataROI[0])
            if event.ydata > max(self.dataROI[1]):
                event.ydata = max(self.dataROI[1]) * 0.99999999999
                # max*0.9999 Because max data not accepted with
                # the fitting method.
            if event.ydata < min(self.dataROI[1]):
                event.ydata = 0.0
                # Beause min data can be slightly under 0 due
                # to baseline correction.

            # Add the created peak values to the end of each lists
            self.height = np.concatenate([self.height, [event.ydata]])
            self.wavenb = np.concatenate([self.wavenb, [event.xdata]])
            self.fwhm = np.concatenate([self.fwhm, [5]])
            self.heightrange = np.concatenate([self.heightrange, [np.Inf]])
            self.wavenbrange = np.concatenate([self.wavenbrange, [np.Inf]])
            self.fwhmrange = np.concatenate([self.fwhmrange, [np.Inf]])

            if self.refreshcount != 0:  # If deconv. already done before
                plt.gca().get_lines().pop(-1).remove()  # Remove fitting result
                self.refreshcount = 0
            self.parent.ax.plot(
                self.dataROI[0],
                self.function(
                    self.dataROI[0], event.ydata, 5. / 2, event.xdata),
                linewidth=1.0, picker=pickersize)  # Plots the new peak
            lines = plt.gca().get_lines()
            [plt.setp(line, 'linewidth', 1.0) for line in lines]
            # Default thickness for each peak

            self.lineid = len(lines) - 2  # Selected peak is the one created
            thisline = lines[self.lineid + 1]
            plt.setp(thisline, 'linewidth', 2.0)  # Created peak thicker

            self.peakcreation = True  # Activation of the peake creation state
            # (deactivates the automatic update of the values due to the entry
            # text change)
            self.heightsel.set(str(event.ydata))
            self.wavenbsel.set(str(event.xdata))
            self.fwhmsel.set(str(5))
            self.heightrangesel.set(str(np.Inf))
            self.wavenbrangesel.set(str(np.Inf))
            self.fwhmrangesel.set(str(np.Inf))
            self.peakcreation = False  # Deactivation of peak creation state

            # Increments the number of peaks
            self.nbpeak += 1
            self.peaklab.config(text='Number of Peaks: {}'.format(self.nbpeak))

            self.parent.canvas.draw()  # Updates canvas

    def on_pick(self, event):
        """
        Curve selecting callback.
        Left click = Selection.
        """
        # If no tools from toolbar selected and click = left click
        if event.mouseevent.button == 1 and self.parent.toolbar._active is None:
            # Every lines get their default thickness except the selected one
            lines = plt.gca().get_lines()
            [plt.setp(line, 'linewidth', 1.0) for line in lines]
            thisline = event.artist
            plt.setp(thisline, 'linewidth', 2.0)

            # Searching the selected line index
            lineid = None
            for i in range(len(lines)):
                if lines[i] is thisline:
                    lineid = i - 1

            self.peakcreation = True  # Activation of the peake creation state
            # (deactivates the automatic update of the values due to the entry
            # text change)
            # Ensures selected curve is a peak
            if lineid is not None and lineid != -1 \
                    and lineid < len(self.height):
                # Updates entries values
                self.lineid = lineid
                self.heightsel.set(str(self.height[lineid]))
                self.wavenbsel.set(str(self.wavenb[lineid]))
                self.fwhmsel.set(str(self.fwhm[lineid]))
                self.heightrangesel.set(str(self.heightrange[lineid]))
                self.wavenbrangesel.set(str(self.wavenbrange[lineid]))
                self.fwhmrangesel.set(str(self.fwhmrange[lineid]))
            else:  # No peak selected
                # Updates entries values
                self.heightsel.set('')
                self.wavenbsel.set('')
                self.fwhmsel.set('')
                self.heightrangesel.set('')
                self.wavenbrangesel.set('')
                self.fwhmrangesel.set('')
            self.peakcreation = False  # Deactivation of peak creation state

        self.parent.canvas.draw()  # Updates canvas

    def _quit(self):
        """ Properly quits the Manual Analysis tool """
        self.parent.wm_title('Spectral Analysis')
        for cid in self.eventId:
            self.parent.canvas.mpl_disconnect(cid)
        self.parent.unbind("<Return>")
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.destroy()
