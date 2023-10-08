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

import tkinter as tk
from tkinter.filedialog import askopenfilename
import spc
import numpy as np
import matplotlib.pyplot as plt
from spload import spload2


def find_nearest(array_xy, value_xy, aspect_ratio):
    """
    Return the index of the nearest value in array_xy from value_xy.
    Designed to graphically take the nearest value
    """
    x = array_xy[0]
    y = array_xy[1]
    x0 = value_xy[0]
    y0 = value_xy[1]
    a = (max(x) - min(x)) / (max(y) - min(y)) * aspect_ratio
    y = a * y
    y0 = a * y0
    idx = np.sqrt((x - x0)**2 + (y - y0)**2).argmin()
    return idx


def removal(spectrumData, remData, peakBounds):
    """
    Removes the remData (spectrum of the component to be removed) component inside spectrumData (spectrum of the
    component composed of both the material of interest and the material to be removed from the spectrum), based on the
    main peak feature (highest peak) located inside the peakBounds of remData.
    :param spectrumData: list of [xData (list), yData (list)] of the spectrum of the material of interest.
    :param remData: list of [xData (list), yData (list)] of the spectrum of the material to be removed from the spectrum
        of the material of interest.
    :param peakBounds: list of [x0, x1] containing the boundaries between which to find the highest peak in the spectrum
        of remData.
    :return: list of [xData (list), yData (list)] containing the spectrum of the material of interest removed from the
        component of remData.
    """
    remData = np.array(remData)
    xsortarg = np.argsort(remData[0])
    remData[1] = remData[1][xsortarg]
    remData[0] = remData[0][xsortarg]

    mmaBounds = sorted(peakBounds)
    # Searching the corresponding ID of the min MMA boundary
    mmaInfId = abs(spectrumData[0] - mmaBounds[0]).argmin()
    # Searching the corresponding ID of the max MMA boundary
    mmaSupId = abs(spectrumData[0] - mmaBounds[1]).argmin()

    # Restraining the analysis area to the MMA peak
    xorig = spectrumData[0][min(mmaInfId, mmaSupId):max(mmaInfId, mmaSupId)]
    yorig = spectrumData[1][min(mmaInfId, mmaSupId):max(mmaInfId, mmaSupId)]

    ymmaorig = np.array(remData[1])

    ymma = np.array(ymmaorig[min(mmaInfId, mmaSupId):max(mmaInfId, mmaSupId)])

    # MMA removal algorithm
    coef = 0.  # Sizing coefficient for the MMA spectrum
    coefsave = coef  # Saved sizing coef
    delta = 0.5  # Initial step
    area = 0  # Area defining how "big" is the presence of the MMA peak
    areasave = np.Inf  # Saved area
    n = 2  # Initial loop's number
    # Terminal precision for coef (max of abs(coef(n-1)-coef(n)) wanted)
    precision = 0.0001
    maxloop = 1000  # Terminal loop if precision is not reached
    while(abs(delta) > precision and n < maxloop):
        # 6th order polynomial fit coefficients of the datas
        # (for a smoothed, perfect derivative)
        p = np.polyfit(xorig, yorig - ymma * coef, 6)
        # 1st derivative coefficients
        p1 = [p[i] * (6 - i) for i in range(7)][0:6]
        # 2nd derivative coefficients
        p2 = [p1[i] * (5 - i) for i in range(6)][0:5]
        yfit2 = np.poly1d(p2)(xorig)  # 2nd derivative

        # Calculating the area of positive sign
        for i in range(len(yfit2)):
            if np.sign(yfit2[i]) == 1:  # Should always be positive
                yfit2[i] = 0
        area = abs(np.trapz(yfit2, xorig))

        # Checking if current coef improves the peak removal
        if areasave >= area:
            coefsave = coef
            areasave = area
        # If no improvements, the step change its sign and is reduced
        else:
            # Previous change is cancelled since there was no improvement
            coef -= delta
            # The step is halved and take the opposite sign
            delta = -delta / 2

        coef += delta  # Coef is incremented of delta
        n += 1  # Number of loop is incremented

    # Applying the MMA peak removal with the best coefficient
    data = [spectrumData[0], spectrumData[1] - ymmaorig * coefsave]
    return data


class spectrumRem(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)  # Frame initialization
        self.parent = parent

        # Asking the spectrum's file to subtract
        self.mmaFilePath = askopenfilename(
            title='Open .spc spectrum file to subtract',
            filetypes=[('SP file', '*.sp'), ('SPC file', '*.spc'), ('SPC file', '*.SPC')])

        # Opening the file
        if self.mmaFilePath[-4:] == '.spc':
            self.mmaSpectrum = spc.File(self.mmaFilePath)
            self.mmaDatas = [np.array(self.mmaDatas[0]), np.array(self.mmaDatas[1])]
        elif self.mmaFilePath[-3:] == '.sp':
            y, x = spload2(self.mmaFilePath)
            self.mmaDatas = [np.array(x), np.array(y)]

        # Creating interface
        lab = tk.Label(master=self,
                       text='Reference peak boundaries (wave numbers)')
        lab.grid(row=0, column=0)
        self.binf_var = tk.StringVar(value=self.parent.bounds_sub[list(self.parent.bounds_sub.keys())[0]][0])
        self.binf_var.trace('w', self.binf_callback)
        self.bsup_var = tk.StringVar(value=self.parent.bounds_sub[list(self.parent.bounds_sub.keys())[0]][1])
        self.bsup_var.trace('w', self.bsup_callback)
        self.binfentry = tk.Entry(master=self, textvariable=self.binf_var)
        self.binfentry.grid(row=0, column=1)
        self.bsupentry = tk.Entry(master=self, textvariable=self.bsup_var)
        self.bsupentry.grid(row=0, column=2)
        self.okbut = tk.Button(master=self, text='Done',
                          command=self.subtract)
        self.okbut.grid(row=0, column=3, rowspan=10, sticky=tk.NE+tk.SW)

        self.absorb_lab = tk.Label(self, text='Transmittance to absorbance')
        self.absorb_lab.grid(row=1, column=0, sticky=tk.E)
        self.absorb_var = tk.IntVar()
        self.absorb_check = tk.Checkbutton(self, variable=self.absorb_var, command=self.absorb)
        self.absorb_check.grid(row=1, column=1)

        self.parent.bind("<Return>", lambda event: self.apply())

        # Displaying the spectrum to subtract
        self.mmaLine = self.parent.ax.plot(self.mmaDatas[0], self.mmaDatas[1], 'k-', label='Spectrum to remove')[0]
        plt.legend()

        # Visual representation of the boundaries
        self.x1 = None  # 1st extremity of the boundaries
        self.x2 = None  # 2nd extremity of the boundaries
        self.nbpts = 0  # Check variable (how many points have been created)
        self.ptnum = 0  # Number of the point selected
        self.click = False  # Check variable (is mouse button pressed)
        self.motion = False


        # Getting the original values based on the config file
        self.x1 = abs(self.parent.data[0] - float(self.binfentry.get())).argmin()
        # Plotting the x1 point
        self.gx1 = self.parent.ax.plot(
            self.parent.data[0][self.x1], self.parent.data[1][self.x1],
            'ro', picker=10)[0]

        self.x2 = abs(self.parent.data[0] - float(self.bsupentry.get())).argmin()
        # Plotting the x2 point
        self.gx2 = self.parent.ax.plot(
            self.parent.data[0][self.x2], self.parent.data[1][self.x2],
            'ro', picker=10)[0]

        # Plotting the x1-x2 line
        self.gx1x2 = self.parent.ax.plot(
            self.parent.data[0][[self.x1, self.x2]],
            self.parent.data[1][[self.x1, self.x2]], 'r-')[0]

        # Connecting the matplotlib events
        self.eventId = []
        self.eventId.append(self.parent.canvas.mpl_connect('button_press_event', self.on_click))
        self.eventId.append(self.parent.canvas.mpl_connect('button_release_event', self.on_release))
        self.eventId.append(self.parent.canvas.mpl_connect('pick_event', self.on_pick))
        self.eventId.append(self.parent.canvas.mpl_connect('motion_notify_event', self.on_motion))

    def absorb(self):
        if self.absorb_var.get():
            self.mmaDatas[1] = -np.log10(self.mmaDatas[1]/100)
        else:
            self.mmaDatas[1] = 100*np.exp(-self.mmaDatas[1]*np.log(10))
        self.mmaLine.set_ydata(self.mmaDatas[1])

        # Change the graphical boundaries of the subtracting peak
        self.gx1.set_data(self.parent.data[0][self.x1],
                          self.parent.data[1][self.x1])
        self.gx2.set_data(self.parent.data[0][self.x2],
                          self.parent.data[1][self.x2])
        if self.x1 is not None and self.x2 is not None:
            self.gx1x2.set_data(
                self.parent.data[0][[self.x1, self.x2]],
                self.parent.data[1][[self.x1, self.x2]])

        # recompute the ax.dataLim
        self.parent.ax.relim()
        # update ax.viewLim using the new dataLim
        self.parent.ax.autoscale_view()
        if self.parent.display:
            self.parent.canvas.draw()

    def subtract(self):
        mmaBoundsInf = float(self.binf_var.get())
        mmaBoundsSup = float(self.bsup_var.get())
        mmaBounds = sorted([mmaBoundsInf, mmaBoundsSup])
        mmaData = [self.mmaDatas[0], self.mmaDatas[1]]
        self.parent.data_change(removal(self.parent.data, mmaData, mmaBounds))
        self._quit()

    def on_motion(self, event):
        """
        Mouse motion callback. Drags the selected point and chenge the
        associated boundary value
        """
        # If mouse pressed and a point is created and no tools from
        # toolbar are selected
        self.motion = True
        if self.click is True and self.ptnum != 0 \
                and self.parent.toolbar._active is None:
            # Aspect ratio for better mouse selection
            aspect_ratio = float(self.parent.winfo_height()) / \
                self.parent.winfo_width()

            # Mouse click coordinates
            x = event.xdata
            y = event.ydata

            # Nearest point on curve
            xid = find_nearest([self.parent.data[0], self.parent.data[1]],
                               [x, y], aspect_ratio)

            # Storing the selected boundary, plotting the point and the
            # linking line
            if self.ptnum == 1:
                self.x1 = xid
                self.gx1.set_data(self.parent.data[0][xid],
                                                  self.parent.data[1][xid])
            else:
                self.x2 = xid
                self.gx2.set_data(self.parent.data[0][xid],
                                                  self.parent.data[1][xid])
            if self.x1 is not None and self.x2 is not None:
                self.gx1x2.set_data(
                    self.parent.data[0][[self.x1, self.x2]],
                    self.parent.data[1][[self.x1, self.x2]])

            # Updates the entry widgets
            self.binf_var.set(str(self.parent.data[0][self.x1]))
            self.bsup_var.set(str(self.parent.data[0][self.x2]))
            # Updates the canvas
            if self.parent.display:
                self.parent.canvas.draw()
            self.motion = False

    def binf_callback(self, *args):
        if self.motion is False:
            x = float(self.binfentry.get())

            # Nearest point on curve
            xid = abs(self.parent.data[0] - x).argmin()

            # Storing the selected boundary, plotting the point and the
            # linking line
            self.x1 = xid
            self.gx1.set_data(self.parent.data[0][xid],
                                                  self.parent.data[1][xid])

            if self.x1 is not None and self.x2 is not None:
                self.gx1x2.set_data(
                    self.parent.data[0][[self.x1, self.x2]],
                    self.parent.data[1][[self.x1, self.x2]])

            # Updates the canvas
            if self.parent.display:
                self.parent.canvas.draw()

    def bsup_callback(self, *args):
        if self.motion is False:
            x = float(self.bsupentry.get())

            # Nearest point on curve
            xid = abs(self.parent.data[0] - x).argmin()

            # Storing the selected boundary, plotting the point and the
            # linking line
            self.x2 = xid
            self.gx2.set_data(self.parent.data[0][xid],
                                                  self.parent.data[1][xid])

            if self.x1 is not None and self.x2 is not None:
                self.gx1x2.set_data(
                    self.parent.data[0][[self.x1, self.x2]],
                    self.parent.data[1][[self.x1, self.x2]])

            # Updates the canvas
            if self.parent.display:
                self.parent.canvas.draw()

    def on_release(self, event):
        """ Mouse button released callback. Deactivates the pressed state """
        self.click = False
        self.ptnum = 0

    def on_click(self, event):
        """
        Activate the button pressed state (self.click = True)
        """
        # Activates the button pressed state
        self.click = True

    def on_pick(self, event):
        """ Curve selecting callback """
        if self.parent.toolbar._active is None:  # If no tool from toolbar is selected

            # Getting the selected line
            lines = [self.gx1, self.gx2]
            thisline = event.artist
            if thisline is self.gx1:
                self.ptnum = 1
            if thisline is self.gx2:  # If second point
                self.ptnum = 2
            # Restores other lines thickness to default
            [plt.setp(line, 'mew', 0.5) for line in lines]
            # Increase thickness of selected line
            plt.setp(thisline, 'mew', 2.0)

            # Update canavas
            if self.parent.display:
                self.parent.canvas.draw()

    def _quit(self):
        """ Properly quits the spectrum removal tool """
        self.parent.wm_title('Spectral Analysis')
        for i in self.eventId:
            self.parent.canvas.mpl_disconnect(i)
        self.parent.unbind("<Return>")
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.destroy()
        return
