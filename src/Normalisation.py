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


class Normalisation(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)  # Frame initialization
        self.parent = parent

        # Visual representation of the reference peak
        self.x1 = None  # 1st extremity of the boundaries
        self.click = False  # Check variable (is mouse button pressed)
        self.motion = False
        self.picked = False

        # Getting the original values based on the max peak
        self.x1 = self.parent.data[1].argmax()
        # Plotting the x1 point
        self.gx1 = self.parent.ax.plot(
            self.parent.data[0][self.x1], self.parent.data[1][self.x1],
            'ro', picker=10)[0]

        # Creating interface
        lab = tk.Label(master=self,
                       text='Reference peak (wave numbers)')
        lab.pack(side=tk.LEFT)
        self.refPeakVar = tk.StringVar(value=str(self.parent.data[0][self.x1]))
        self.refPeakVar.trace('w', self.refPeak_callback)
        self.refPeak_entry = tk.Entry(self, textvariable=self.refPeakVar)
        self.refPeak_entry.pack(side=tk.LEFT)
        self.apply_button = tk.Button(self, text='Apply', command=self.apply)
        self.apply_button.pack(side=tk.LEFT)
        self.parent.bind("<Return>", lambda event: self.apply())

        # Connecting the matplotlib events
        self.eventId = []
        self.eventId.append(self.parent.canvas.mpl_connect('button_press_event', self.on_click))
        self.eventId.append(self.parent.canvas.mpl_connect('button_release_event', self.on_release))
        self.eventId.append(self.parent.canvas.mpl_connect('pick_event', self.on_pick))
        self.eventId.append(self.parent.canvas.mpl_connect('motion_notify_event', self.on_motion))

    def on_motion(self, event):
        """
        Mouse motion callback. Drags the selected point and chenge the
        associated boundary value
        """
        # If mouse pressed and a point is created and no tools from
        # toolbar are selected
        self.motion = True
        if self.click is True and self.parent.toolbar._active is None and self.picked is True:
            # Aspect ratio for better mouse selection
            aspect_ratio = float(self.parent.winfo_height()) / \
                self.parent.winfo_width()

            # Mouse click coordinates
            x = event.xdata
            y = event.ydata

            # Nearest point on curve
            xid = find_nearest([self.parent.data[0], self.parent.data[1]],
                               [x, y], aspect_ratio)

            # Storing the selected boundary, plotting the point and the
            # linking line
            self.x1 = xid
            self.gx1.set_data(self.parent.data[0][xid], self.parent.data[1][xid])

            # Updates the entry widgets
            self.refPeakVar.set(str(self.parent.data[0][self.x1]))
            # Updates the canvas
            if self.parent.display:
                self.parent.canvas.draw()
            self.motion = False

    def refPeak_callback(self, *args):
        if self.motion is False:
            x = float(self.bsupentry.get())

            # Nearest point on curve
            xid = abs(self.parent.data[0] - x).argmin()

            # Storing the selected boundary, plotting the point and the
            # linking line
            self.x1 = xid
            self.gx1.set_data(self.parent.data[0][xid], self.parent.data[1][xid])

            # Updates the canvas
            if self.parent.display:
                self.parent.canvas.draw()

    def on_release(self, event):
        """ Mouse button released callback. Deactivates the pressed state """
        self.click = False
        self.picked = False

    def on_click(self, event):
        """
        Activate the button pressed state (self.click = True)
        """
        # Activates the button pressed state
        self.click = True

    def on_pick(self, event):
        """ Curve selecting callback """
        if self.parent.toolbar._active is None:  # If no tool from toolbar is selected
            self.picked = True

    def apply(self):
        trange1p = (max(self.parent.data[0]) - min(self.parent.data[0]))*0.01
        minb = np.abs(self.parent.data[0] - (self.parent.data[0][self.x1]-trange1p)).argmin()
        maxb = np.abs(self.parent.data[0] - (self.parent.data[0][self.x1]+trange1p)).argmin()
        norm_peak = max(self.parent.data[1][min(minb, maxb):max(minb, maxb)])
        self.parent.data_change([self.parent.data[0], self.parent.data[1] / norm_peak])
        self._quit()

    def _quit(self):
        """ Properly quits the spectrum removal tool """
        self.parent.wm_title('Spectral Analysis')
        for i in self.eventId:
            self.parent.canvas.mpl_disconnect(i)
        self.parent.unbind("<Return>")
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.destroy()
        return
