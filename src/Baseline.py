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
from spectrumFit import baseline


class Baseline(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)  # Frame initialization
        self.parent = parent

        # Creating interface
        tk.Label(self, text='Select the order of the baseline : ').pack(side=tk.LEFT)
        self.order_list = tk.Listbox(self, selectmode=tk.EXTENDED)
        self.order_list.pack(side=tk.LEFT)
        self.order_list.insert(tk.END, 'Multi-Linear')
        self.order_list.insert(tk.END, 'Quadratic')
        self.order_list.insert(tk.END, 'Cubic')
        self.order_list.bind('<<ListboxSelect>>', self.select)

        self.apply_button = tk.Button(self, text='Apply', command=self.apply)
        self.apply_button.pack(side=tk.LEFT)

        self.parent.bind("<Return>", lambda event: self.apply())

        # Variable initialization
        self.bl = None
        self.blLine = None

    def select(self, event, *args, **kwargs):
        method = event.widget.get(int(event.widget.curselection()[0]))
        if method == 'Multi-Linear':
            self.bl = baseline(self.parent.data[0], self.parent.data[1], method='lin', ax = self.parent.ax)
        elif method == 'Quadratic':
            self.bl = baseline(self.parent.data[0], self.parent.data[1], method='quad', ax = self.parent.ax)
        elif method == 'Cubic':
            self.bl = baseline(self.parent.data[0], self.parent.data[1], method='cub', ax = self.parent.ax)
        else:
            method = None

        if method is not None:
            if self.blLine is None:
                self.blLine = self.parent.ax.plot(self.parent.data[0], self.bl, 'k', label='Baseline')[0]
                plt.legend()
            else:
                self.blLine.set_ydata(self.bl)
            if self.parent.display:
                self.parent.canvas.draw()

    def apply(self):
        if self.bl is not None:
            self.parent.data_change([self.parent.data[0], self.parent.data[1]-self.bl])
        self._quit()

    def _quit(self):
        """ Properly quits the spectrum removal tool """
        self.parent.wm_title('Spectral Analysis')
        self.parent.unbind("<Return>")
        self.parent.data_change(self.parent.data, refresh_ax_only=True)
        self.destroy()
        return
