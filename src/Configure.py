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
Configure class.
Defines a configuration window to change the various preset parameters for the
automatic analysis.
"""

import tkinter as tk
from functools import partial
from functools import partial


FILENAME = './configure.ini'
RESOLUTION = '648x486'


def paramget(array, paramname):
    """
    Look for the value of the parameter 'paramname' in array.
    First column of array must be the parameter's name.
    """
    if not array:
        return None
    for i in range(len(array)):
        if array[i][0] == paramname:
            return array[i][1:len(array[i])]
    return None


class Configure(tk.Toplevel):

    def __init__(self, parent, *args, **kwargs):
        """ Initialize every parameters for the configuration window """
        self.parent = parent
        tk.Toplevel.__init__(self, parent, *args, **kwargs)  # Configuration window

        self.canvconfwin = tk.Canvas(master=self)
        self.canvconfwin.grid(column=0, row=0, sticky=tk.N+tk.S+tk.E+tk.W)

        self.scrollbarh = tk.Scrollbar(self, command=self.canvconfwin.xview, orient=tk.HORIZONTAL)
        self.scrollbarh.grid(column=0, row=1, sticky=tk.E+tk.W)
        self.canvconfwin.configure(xscrollcommand=self.scrollbarh.set)
        self.scrollbarv = tk.Scrollbar(self, command=self.canvconfwin.yview, orient=tk.VERTICAL)
        self.scrollbarv.grid(column=1, row=0, sticky=tk.N+tk.S)
        self.canvconfwin.configure(yscrollcommand=self.scrollbarv.set)

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.confwin = tk.Frame(self.canvconfwin)
        self.window = self.canvconfwin.create_window((0,0), window=self.confwin, anchor=tk.NW)

        self.canvconfwin.bind_all('<Configure>', self.on_configure)
        self.bind_all('<MouseWheel>', self.on_mousewheel)
        self.bind_all('<Button-4>', self.on_mousewheel)
        self.bind_all('<Button-5>', self.on_mousewheel)


        # Getting infos from the configure.ini file
        self.bounds_sub = self.parent.bounds_sub
        self.areabounds = self.parent.areaBounds
        self.guess = self.parent.guess
        self.bounds = self.parent.peakBounds
        self.orderArea = []  # Keep track of the order of the area
        self.orderAreaSub = []  # Keep track of the order of the area sub

        # GUI creation for modifying the parameters
        row = 0
        col = 0

        # Subtraction boundaries
        self.bounds_sub_labels = {}
        self.bounds_sub_edit = {}
        for i in range(len(self.bounds_sub)):
            self.orderAreaSub.append(self.bounds_sub.keys()[i])
            self.bounds_sub_labels[self.bounds_sub.keys()[i]] = tk.Label(master=self.confwin, text=self.bounds_sub.keys()[i]+' subtract boundaries')
            self.bounds_sub_labels[self.bounds_sub.keys()[i]].grid(row=row, column=0, padx=2, pady=2)

            self.bounds_sub_edit[self.bounds_sub.keys()[i]] = [tk.Entry(master=self.confwin, width=5) for j in range(len(self.bounds_sub[self.bounds_sub.keys()[i]]))]
            [self.bounds_sub_edit[self.bounds_sub.keys()[i]][j].insert(0, str(self.bounds_sub[self.bounds_sub.keys()[i]][j])) for j in range(len(self.bounds_sub[self.bounds_sub.keys()[i]]))]
            [self.bounds_sub_edit[self.bounds_sub.keys()[i]][j].grid(row=row, column=j+1, padx=1, pady=2) for j in range(len(self.bounds_sub[self.bounds_sub.keys()[i]]))]

            row += 1

        # Areas
        self.areabounds_labels = {}
        self.areabounds_edit = {}
        self.guess_labels = {}
        self.guess_edit = {}
        self.boundsinf_labels = {}
        self.boundsinf_edit = {}
        self.boundssup_labels = {}
        self.boundssup_edit = {}
        self.add_buttons = {}
        self.del_buttons = {}
        self.indics_label = {}
        self.rows = {}
        self.keys = []
        l = 0  # Indice corresponding to the area being treated
        padxpeak = [(1,1), (1,1), (1,6)]
        for areaname in self.areabounds:
            self.orderArea.append(areaname)

            # Area boundaries
            self.rows[areaname] = row
            self.areabounds_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' area boundaries')
            self.areabounds_labels[areaname].grid(row=row, column=0, padx=2, pady=(10,1), sticky=tk.E)
            self.areabounds_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.areabounds[areaname]))]
            [self.areabounds_edit[areaname][k].insert(0, str(self.areabounds[areaname][k])) for k in range(len(self.areabounds[areaname]))]
            [self.areabounds_edit[areaname][k].grid(row=row, column=k+1, padx=1, pady=(10,1)) for k in range(len(self.areabounds[areaname]))]
            row += 1

            # Peaks indications
            indics = ['Height', 'Half FWHM', 'Wave Nb']
            self.indics_label[areaname] = [tk.Label(master=self.confwin, text=indics[k%3]+str(k/3), font=("Helvetica", 5)) for k in range(len(self.guess[areaname]))]
            [self.indics_label[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=(1,0)) for k in range(len(self.indics_label[areaname]))]
            row += 1

            # Guess
            self.guess_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' initial guess')
            self.guess_labels[areaname].grid(row=row, column=0, padx=2, pady=1, sticky=tk.E)
            self.guess_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.guess[areaname]))]
            [self.guess_edit[areaname][k].insert(0, str(self.guess[areaname][k])) for k in range(len(self.guess[areaname]))]
            [self.guess_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=1) for k in range(len(self.guess[areaname]))]

            # Add and del peak buttons
            self.add_buttons[areaname] = tk.Button(master=self.confwin, text='Add\npeak', command=partial(self.addpeak, areaname))
            self.add_buttons[areaname].grid(row=row, column=len(self.guess[areaname])+1, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
            self.del_buttons[areaname] = tk.Button(master=self.confwin, text='Del\npeak', command=partial(self.delpeak, areaname))
            self.del_buttons[areaname].grid(row=row, column=len(self.guess[areaname])+3, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
            row += 1

            # Bounds inf
            self.boundsinf_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' inferior boundaries')
            self.boundsinf_labels[areaname].grid(row=row, column=0, padx=2, pady=1, sticky=tk.E)
            self.boundsinf_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.bounds[areaname][0]))]
            [self.boundsinf_edit[areaname][k].insert(0, str(self.bounds[areaname][0][k])) for k in range(len(self.bounds[areaname][0]))]
            [self.boundsinf_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=1) for k in range(len(self.bounds[areaname][0]))]
            row += 1

            # Bounds sup
            self.boundssup_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' superior boundaries')
            self.boundssup_labels[areaname].grid(row=row, column=0, padx=2, pady=(1,3), sticky=tk.E)
            self.boundssup_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.bounds[areaname][1]))]
            [self.boundssup_edit[areaname][k].insert(0, str(self.bounds[areaname][1][k])) for k in range(len(self.bounds[areaname][1]))]
            [self.boundssup_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=(1,3)) for k in range(len(self.bounds[areaname][1]))]
            row += 1

            l += 1
            self.keys.append(areaname)

        self.finishRow = row
        self.nbArea = l-1

        self.addarea_button = tk.Button(master=self.confwin, text='Add Area', command=self.addarea)
        self.addarea_button.grid(row=199, column=1, columnspan=3, sticky=tk.W+tk.E)

        self.delarea_button = tk.Button(master=self.confwin, text='Delete Previous Area', command=self.delarea)
        self.delarea_button.grid(row=199, column=4, columnspan=3, sticky=tk.W+tk.E)

        self.apply_button = tk.Button(master=self.confwin, text='Apply', command=self.apply)
        self.apply_button.grid(row=200, column=1, pady=10, padx=1, columnspan=2)


    def on_configure(self, event):
        self.canvconfwin.configure(scrollregion=self.canvconfwin.bbox('all'))

    def on_mousewheel(self, event):
        if event.num == 4 or event.delta == 120:
            self.canvconfwin.yview('scroll', -1, 'units')
        if event.num == 5 or event.delta == -120:
            self.canvconfwin.yview('scroll', 1, 'units')


    def apply(self):
        # Changing the values in self.conf
        for i in range(len(self.conf)):
            if self.conf[i][0] == 'mmapeakbound':
                self.conf[i][1] = self.mmapeakinfentry.get()
                self.conf[i][2] = self.mmapeaksupentry.get()

        # Writing self.conf as the new configuration file
        with open(FILENAME, 'w+') as conffile:
            for line in self.conf:
                totline = ''
                for col in line:
                    totline += col + ' '
                conffile.write(totline + '\n')

    def addarea(self):
        global newareaname
        askareaname = tk.Toplevel(master=self.confwin)
        areaname_lab = tk.Label(master=askareaname, text='New area name : ')
        areaname_lab.pack(side=tk.LEFT)
        areaname_entry = tk.Entry(master=askareaname)
        areaname_entry.pack(side=tk.LEFT)

        def askareanameOK():
            global newareaname
            newareaname = areaname_entry.get()
            if newareaname == '':
                newareaname = 'New area'
                ok = False
                i = 0
                while ok is False:
                    ok = True
                    i += 1
                    for areaname in self.keys:
                        if areaname == newareaname:
                            ok = False
                            if areaname[-1] == str(i-1):
                                newareaname = newareaname[:-1] + str(i)
                            else:
                                newareaname = newareaname + str(i)
                            break

            askareaname.destroy()

        ok_button = tk.Button(master=askareaname, text='OK', command=askareanameOK)
        ok_button.pack(side=tk.LEFT)

        self.confwin.wait_window(askareaname)

        self.orderArea.append(newareaname)
        self.areabounds[newareaname] = [0.0, 1.0]

        self.rows[newareaname] = self.finishRow
        self.areabounds_labels[newareaname] = tk.Label(master=self.confwin, text=newareaname+' area boundaries')
        self.areabounds_labels[newareaname].grid(row=self.finishRow, column=0, padx=2, pady=(10,1), sticky=tk.E)
        self.areabounds_edit[newareaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.areabounds[newareaname]))]
        [self.areabounds_edit[newareaname][k].insert(0, str(self.areabounds[newareaname][k])) for k in range(len(self.areabounds[newareaname]))]
        [self.areabounds_edit[newareaname][k].grid(row=self.finishRow, column=k+1, padx=1, pady=(10,1)) for k in range(len(self.areabounds[newareaname]))]
        self.finishRow += 1

        areaname = newareaname
        row = self.finishRow
        padxpeak = [(1,1), (1,1), (1,6)]
        indics = ['Height', 'Half FWHM', 'Wave Nb']
        l = self.nbArea + 1

        self.guess[areaname] = [0.0001, 5.0, 0.0001]
        self.indics_label[areaname] = [tk.Label(master=self.confwin, text=indics[k%3]+str(k/3), font=("Helvetica", 5)) for k in range(len(self.guess[areaname]))]
        [self.indics_label[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=(1,0)) for k in range(len(self.indics_label[areaname]))]
        row += 1
        # Guess
        self.guess_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' initial guess')
        self.guess_labels[areaname].grid(row=row, column=0, padx=2, pady=1, sticky=tk.E)
        self.guess_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.guess[areaname]))]
        [self.guess_edit[areaname][k].insert(0, str(self.guess[areaname][k])) for k in range(len(self.guess[areaname]))]
        [self.guess_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=1) for k in range(len(self.guess[areaname]))]

        # Add and del peak buttons
        self.add_buttons[areaname] = tk.Button(master=self.confwin, text='Add\npeak', command=partial(self.addpeak, areaname))
        self.add_buttons[areaname].grid(row=row, column=len(self.guess[areaname])+1, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
        self.del_buttons[areaname] = tk.Button(master=self.confwin, text='Del\npeak', command=partial(self.delpeak, areaname))
        self.del_buttons[areaname].grid(row=row, column=len(self.guess[areaname])+3, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
        row += 1

        # Bounds inf
        self.bounds[areaname] = [[0.0001, 0.0001, 0.0001], [1.0, 20.0, 1000.0]]
        self.boundsinf_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' inferior boundaries')
        self.boundsinf_labels[areaname].grid(row=row, column=0, padx=2, pady=1, sticky=tk.E)
        self.boundsinf_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.bounds[areaname][0]))]
        [self.boundsinf_edit[areaname][k].insert(0, str(self.bounds[areaname][0][k])) for k in range(len(self.bounds[areaname][0]))]
        [self.boundsinf_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=1) for k in range(len(self.bounds[areaname][0]))]
        row += 1

        # Bounds sup
        self.boundssup_labels[areaname] = tk.Label(master=self.confwin, text=areaname+' superior boundaries')
        self.boundssup_labels[areaname].grid(row=row, column=0, padx=2, pady=(1,3), sticky=tk.E)
        self.boundssup_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(len(self.bounds[areaname][1]))]
        [self.boundssup_edit[areaname][k].insert(0, str(self.bounds[areaname][1][k])) for k in range(len(self.bounds[areaname][1]))]
        [self.boundssup_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=(1,3)) for k in range(len(self.bounds[areaname][1]))]
        row += 1

        l += 1

        self.finishRow = row
        self.nbArea = l - 1
        self.keys.append(areaname)
        return

    def addpeak(self, areaname, *args):
        #self.guess[areaname] = np.concatenate() append([0.0001, 5.0, 0.0001])

        row = self.rows[areaname]+1
        padxpeak = [(1,1), (1,1), (1,6)]
        indics = ['Height', 'Half FWHM', 'Wave Nb']
        maxcol = len(self.indics_label[areaname]) -1
        self.indics_label[areaname] += [tk.Label(master=self.confwin, text=indics[k%3]+str(k/3), font=("Helvetica", 5)) for k in range(maxcol+1, maxcol+4)]
        [self.indics_label[areaname][k].grid(row=row, column=k, padx=padxpeak[k%3], pady=(1,0)) for k in range(maxcol+1, maxcol+4)]

        # Guess
        self.guess_edit[areaname] += [tk.Entry(master=self.confwin, width=5) for k in range(maxcol+1, maxcol+4)]
        [self.guess_edit[areaname][k].insert(0, str(self.guess[areaname][k])) for k in range(maxcol+1, maxcol+4)]
        [self.guess_edit[areaname][k].grid(row=row, column=k, padx=padxpeak[k%3], pady=1) for k in range(maxcol+1, maxcol+4)]

        # Add and del peak buttons
        self.add_buttons[areaname].grid(row=row, column=maxcol+4, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
        self.del_buttons[areaname].grid(row=row, column=maxcol+6, padx=1, pady=1, rowspan=3, columnspan=2, sticky=tk.N+tk.S+tk.W)
        row += 1

        # Bounds inf
        self.bounds[areaname] = [[0.0001, 0.0001, 0.0001], [1.0, 20.0, 1000.0]]
        self.boundsinf_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(maxcol+1, maxcol+4)]
        [self.boundsinf_edit[areaname][k].insert(0, str(self.bounds[areaname][0][k])) for k in range(maxcol+1, maxcol+4)]
        [self.boundsinf_edit[areaname][k].grid(row=row, column=k, padx=padxpeak[k%3], pady=1) for k in range(maxcol+1, maxcol+4)]
        row += 1

        # Bounds sup
        self.boundssup_edit[areaname] = [tk.Entry(master=self.confwin, width=5) for k in range(maxcol+1, maxcol+4)]
        [self.boundssup_edit[areaname][k].insert(0, str(self.bounds[areaname][1][k])) for k in range(maxcol+1, maxcol+4)]
        [self.boundssup_edit[areaname][k].grid(row=row, column=k+1, padx=padxpeak[k%3], pady=(1,3)) for k in range(maxcol+1, maxcol+4)]
        row += 1

        l += 1

        self.finishRow = row
        return

    def delarea(self):
        areaname = self.orderArea[-1]
        print(areaname)
        self.boundssup_labels[areaname].destroy()
        [self.boundssup_edit[areaname][i].destroy() for i in range(len(self.boundssup_edit[areaname]))]
        self.bounds.pop(areaname)
        self.boundssup_edit.pop(areaname)
        self.boundssup_labels.pop(areaname)
        self.boundsinf_labels[areaname].destroy()
        [self.boundsinf_edit[areaname][i].destroy() for i in range(len(self.boundsinf_edit[areaname]))]
        self.boundsinf_labels.pop(areaname)
        self.boundsinf_edit.pop(areaname)
        self.add_buttons[areaname].destroy()
        self.add_buttons.pop(areaname)
        self.del_buttons[areaname].destroy()
        self.del_buttons.pop(areaname)
        self.guess_labels[areaname].destroy()
        [self.guess_edit[areaname][i].destroy() for i in range(len(self.guess_edit[areaname]))]
        self.guess_labels.pop(areaname)
        self.guess_edit.pop(areaname)
        self.guess.pop(areaname)
        self.areabounds.pop(areaname)
        self.areabounds_labels[areaname].destroy()
        [self.areabounds_edit[areaname][i].destroy() for i in range(len(self.areabounds_edit[areaname]))]
        self.areabounds_labels.pop(areaname)
        self.areabounds_edit.pop(areaname)
        [self.indics_label[areaname][i].destroy() for i in range(len(self.indics_label[areaname]))]
        self.indics_label.pop(areaname)
        self.rows.pop(areaname)
        self.finishRow -= 5
        self.nbArea -= 1
        self.keys.pop()
        self.orderArea.pop()
        return

    def delpeak(self, areaname, *args):
        return
