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
ROI selection class
Defines a ROI selection window.
Keep the boundaries in instance variables to be used for other functionalities.
"""

import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import spectrumFit as sf
from copy import deepcopy

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


class ROIselect(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        """
        Creates the specific windows in the main parent window and takes
        the figure, canvas, axe and toolbar to be used.
        """
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.parent.wm_title('Spectral Analysis | ROI SELECTION')
        self.x1 = None  # 1st extremity of the boundaries
        self.x2 = None  # 2nd extremity of the boundaries
        self.nbpts = 0  # Check variable (how many points have been created)
        self.ptnum = 0  # Number of the point selected
        self.click = False  # Check variable (is mouse button pressed)

        # Connecting the matplotlib events
        self.eventId = []
        self.eventId.append(self.parent.canvas.mpl_connect('button_press_event', self.on_click))
        self.eventId.append(self.parent.canvas.mpl_connect('button_release_event', self.on_release))
        self.eventId.append(self.parent.canvas.mpl_connect('pick_event', self.on_pick))
        self.eventId.append(self.parent.canvas.mpl_connect('motion_notify_event', self.on_motion))

        # ROI description
        self.roi_desc = tk.Label(self, text="Click on the curve to define your two boundaries.\n(You can move them afterward)")
        self.roi_desc.pack(side=tk.TOP)

        # ROI applying
        self.ok_butt = tk.Button(self, text='Apply', command=self.apply)
        self.ok_butt.pack(side=tk.TOP)

        self.parent.bind("<Return>", lambda event: self.apply())

    def on_click(self, event):
        """
        Mouse click callback. Creates point symbolising a boundary (no more
        than 2 points created). Activate the button pressed state
        (self.click = True)
        """
        # If no more than 2 points are already created and if no tools from
        # toolbar is selected
        if self.nbpts < 2 and self.parent.toolbar._active is None:
            aspect_ratio = float(self.parent.winfo_height()) / \
                self.parent.winfo_width()
            # Aspect ratio for better mouse selection

            # Mouse click coordinates
            x = event.xdata
            y = event.ydata

            # Nearest point on curve
            xid = find_nearest([self.parent.data[0], self.parent.data[1]],
                               [x, y], aspect_ratio)

            # Stocking the boundaries
            if self.x1 is None:
                self.x1 = xid
                # Plotting the points
                self.gx1 = self.parent.ax.plot(
                    self.parent.data[0][xid], self.parent.data[1][xid],
                    'ro', picker=10)[0]
                self.ptnum = 1
            else:
                self.x2 = xid
                # Plotting the points
                self.gx2 = self.parent.ax.plot(
                    self.parent.data[0][xid], self.parent.data[1][xid],
                    'ro', picker=10)[0]
                self.ptnum = 2

            # If the 2 points are created, plots the line between them
            if self.x2 is not None:
                self.gx1x2 = self.parent.ax.plot(
                    self.parent.data[0][[self.x1, self.x2]],
                    self.parent.data[1][[self.x1, self.x2]], 'r-')[0]

            self.nbpts += 1  # Number of points incrementation
            self.parent.canvas.draw()  # Canvas updating

        # Activates the button pressed state
        self.click = True

    def on_motion(self, event):
        """
        Mouse motion callback. Drags the selected point and chenge the
        associated boundary value
        """
        # If mouse pressed and a point is created and no tools from
        # toolbar are selected
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

            # Storing the selected boundary, plotting the point and the
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

            # Updates the canvas
            self.parent.canvas.draw()

    def on_release(self, event):
        """ Mouse button released callback. Deactivates the pressed state """
        self.click = False
        self.ptnum = 0

    def on_pick(self, event):
        """ Curve selecting callback """
        if self.parent.toolbar._active is None:  # If no tool from toolbar is selected

            # Getting the selected line
            lines = [self.gx1, self.gx2]
            thisline = event.artist
            if thisline is self.gx1:
                self.ptnum = 1
            if thisline is self.gx2:  # If second point
                self.ptnum = 2
            # Restores other lines thickness to default
            [plt.setp(line, 'mew', 0.5) for line in lines]
            # Increase thickness of selected line
            plt.setp(thisline, 'mew', 2.0)

            # Update canavas
            self.parent.canvas.draw()

    def apply(self):
        """ Modify the datas to match the selected ROI """
        if self.x1 is not None and self.x2 is not None:
            roiId = sorted([self.x1, self.x2])
            self.parent.bounds = roiId
            data = deepcopy(self.parent.data)
            data[1][roiId[0]:roiId[1]] -= sf.baseline(data[0][roiId[0]:roiId[1]], data[1][roiId[0]:roiId[1]])
            self.parent.data_change(data)
            self._quit()
            self.parent.curWin = None

    def _quit(self):
        """ Properly quits the ROI selection mode """
        self.parent.wm_title('Spectral Analysis')
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.parent.unbind("<Return>")
        self.destroy()
        return
