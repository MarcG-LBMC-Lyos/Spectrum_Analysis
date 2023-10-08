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

import numpy as np
import scipy.optimize as opt
import random as rand
from scipy.special import wofz
from scipy.stats import linregress
import warnings
warnings.simplefilter('ignore', np.RankWarning)


def lorentzianSum(w, *args):
    """ Lorentzian defined by f(w) = L0*wL**2/((w-w0)**2 + wL**2)
    L0 : Max height of the curve
    wL : Half FWHM
    w0 : Position of the max height
    args[i % 3] = L0, args[i % 3 + 1] = FWHM/2 = wL, args[i % 3 + 2] = w0 """
    s = 0
    for i in range(0, len(args), 3):
        s += args[i]*args[i+1]**2/((w-args[i+2])**2+args[i+1]**2)
    return s


def gaussianSum(v, *args):
    """ Gaussian defined by f(v) = G0*exp(-((v-v0)/wG)**2)
    G0 : Max height of the curve
    wG : FWHM/(2*sqrt(2))
    v0 : Position of the max height
    args[i] = G0, args[i+1] = FWHM/2, args[i+2] = v0 """
    s = 0
    for i in range(0, len(args), 3):
        fG = args[i+1]/np.sqrt(np.log(2))
        s += args[i]*np.exp(-((v-args[i+2])/fG)**2)
    return s


def voigtSum(v, *args):
    """ Voigt defined by f(v) = G0*Re(w(z))) with Re = real part function,
    w = Faddeeva function, and z = (v-v0 + i*wL)/wG
    G0 : Max height of the curve
    wG : Half FWHM of Gaussian
    v0 : Position of the max height
    wL : Half FWHM of Lorentzian (= wG here)
    args[i] = G0, args[i+1] = wG, args[i+2] = v0 """
    s = 0
    for i in range(0, len(args), 3):
        wG = args[i+1]
        stemp = np.real(wofz((v-args[i+2] + 1j*wG)))
        s += args[i]*stemp/max(stemp)
    return s


def spectrumFit(f, x, y, **kwargs):
    """ Deconvolve a spectrum by fitting th function f to the spectrum defined
    by (x,y) with optional restriction :
        - guess : the original value of the parameters (write "None" for the
        undefined parameters). The position of the peak will be searched at the
        value +/- 20% of the total length of the x-axis.
        - bounds : for the boundary of the parameters
        - nbPeak : for the total number of peaks in the spectrum """
    rand.seed()
    guess_default = np.array([max(y)/2, (max(x)-min(x))/4, (max(x)+min(x))/2])
    bounds_default_inf = np.array([0, 0, min(x)])
    maxcurve = kwargs.get('maxcurve', max(y))
    bounds_default_sup = np.array([maxcurve, (max(x)-min(x))/2, max(x)])
    guess = kwargs.get('guess', None)
    bounds = kwargs.get('bounds', None)
    nbPeak = kwargs.get('nbPeak', None)
    window = kwargs.get('window', None)

    if window is None:
        window = np.Inf

    if window is not None and (type(window) is not list and type(window) is not np.ndarray):
        window = [window, window, window]

    # If only the position of the peaks are given in "guess"
    if guess is not None:
        guess = np.array(guess)
        guessnoneid = np.argwhere(guess is None)
        for i in range(len(guess)):
            if guess[i] is None:
                guess[i] = bounds_default_inf[i % 3] + \
                    (bounds_default_sup[i % 3]-bounds_default_inf[i % 3])*rand.random()

    if bounds is not None:
        for i in range(len(bounds[0])):
            if bounds[0][i] is None:
                bounds[0][i] = bounds_default_inf[i % 3]
            if bounds[1][i] is None:
                bounds[1][i] = bounds_default_sup[i % 3]

    # User doesn't give any prediction
    if guess is None and bounds is None and nbPeak is None:
        bounds = (bounds_default_inf, bounds_default_sup)
        guess = guess_default
        minerr = np.inf
        poptsave = []
        for i in range(1, 15):
            popt, pcov = opt.curve_fit(f, x, y, p0=guess, bounds=bounds)
            yfit = f(x, *popt)
            erravg = sum(abs(yfit-y))/len(y)
            if erravg < minerr:
                minerr = erravg
                poptsave = popt
            bounds = (
                np.concatenate([bounds[0], bounds_default_inf]),
                np.concatenate([bounds[1], bounds_default_sup]))
            guess = np.concatenate([guess, guess_default])
        return poptsave

    # User gives only the number of peaks in the spectrum
    if guess is None and bounds is None and nbPeak is not None:
        bounds = (bounds_default_inf, bounds_default_sup)
        randnb = rand.random()
        guess = np.array([bounds_default_inf[0] + (bounds_default_sup[0]-bounds_default_inf[0])*randnb,
                 bounds_default_inf[1] + (bounds_default_sup[1]-bounds_default_inf[1])*randnb,
                 bounds_default_inf[2] + (bounds_default_sup[2]-bounds_default_inf[2])*randnb])
        for i in range(1, nbPeak):
            bounds = (
                np.concatenate([bounds[0], bounds_default_inf]),
                np.concatenate([bounds[1], bounds_default_sup]))
            randnb = rand.random()
            guess = np.concatenate([guess,
                np.array([bounds_default_inf[0] + (bounds_default_sup[0]-bounds_default_inf[0])*randnb,
                 bounds_default_inf[1] + (bounds_default_sup[1]-bounds_default_inf[1])*randnb,
                 bounds_default_inf[2] + (bounds_default_sup[2]-bounds_default_inf[2])*randnb])])
        poptsave, pcov = opt.curve_fit(f, x, y, p0=guess, bounds=bounds)
        return poptsave

    # User gives an approximate value of the peaks parameters
    # (mainly the position)
    if guess is not None and bounds is None:
        guess = np.array(guess)
        bounds = (bounds_default_inf, bounds_default_sup)
        for i in range(1, len(guess)/3):
            bounds = (
                np.concatenate([bounds[0], bounds_default_inf]),
                np.concatenate([bounds[1], bounds_default_sup]))
        for i in range(0, len(bounds[0]), 3):
            for j in [0, 1, 2]:
                breakcond = False
                for k in guessnoneid:
                    if i+j == k:
                        breakcond = True
                        break
                if breakcond is True:
                    continue
                bounds[0][i+j] = guess[i+j]*(1-window[j])
                bounds[1][i+j] = guess[i+j]*(1+window[j])
                if bounds[0][i+j] < bounds_default_inf[j]:
                    bounds[0][i+j] = bounds_default_inf[j]
                if bounds[1][i+j] > bounds_default_sup[j]:
                    bounds[1][i+j] = bounds_default_sup[j]
        poptsave, pcov = opt.curve_fit(f, x, y, p0=guess, bounds=bounds)
        return poptsave

    # User gives strict bounds for the peaks
    if bounds is not None:
        bounds = (np.array(bounds[0]), np.array(bounds[1]))
        for i in range(0, len(bounds[0]), 3):
            for j in [0, 1, 2]:
                if bounds[0][i+j] < bounds_default_inf[j]:
                    bounds[0][i+j] = bounds_default_inf[j]
                if bounds[1][i+j] > bounds_default_sup[j]:
                    bounds[1][i+j] = bounds_default_sup[j]
        if guess is None:
            guess = (bounds[0] + bounds[1])/2
        poptsave, pcov = opt.curve_fit(f, x, y, p0=guess, bounds=bounds)
        return poptsave


def derivate(x, y):
    """ Differentiate the curve y = f(x) """
    deriv = np.zeros(len(x))
    deriv[0] = (y[1]-y[0])/(x[1]-x[0])
    deriv[len(x)-1] = (y[len(x)-1]-y[len(x)-2])/(x[len(x)-1]-x[len(x)-2])
    for i in range(1, len(x)-1):
        deriv[i] = \
            ((y[i]-y[i-1])/(x[i]-x[i-1]) + (y[i+1]-y[i])/(x[i+1]-x[i]))/2
    return deriv


def baseline(x, y, **kwargs):
    """ Applies a quadratic baseline correction on a curve """
    method = kwargs.get('method', 'lin')

    idx = locMin(x, y)
    ax = kwargs.get('ax', None)
    if ax is not None:
        ax.plot(x[idx], y[idx], 'or')

    if(method == 'quad'):
        zsave = [0, 0, 0]
        negcountsave = np.inf
        for i in range(0, len(idx)):
            for ii in [0, idx[0]]:
                for iii in [idx[len(idx)-1], len(x)-1]:
                        z = np.polyfit(
                            x[[ii, idx[i], iii]], y[[ii, idx[i], iii]], 2)
                        p = np.poly1d(z)
                        negcount = 0
                        for j in idx:
                            if y[j] - p(x)[j] > 0:
                                negcount += abs(y[j] - p(x)[j])
                            if y[j] - p(x)[j] < 0:
                                if abs(y[j] - p(x)[j]) > 0.01*(max(y)-min(y)):
                                    negcount += np.inf
                                negcount += 40*abs(y[j] - p(x)[j])
                        if z[0] >= 0 and negcount < negcountsave:
                            negcountsave = negcount
                            zsave = z

        baseline = np.poly1d(zsave)(x)
        return baseline

    if(method == 'cub'):
        zsave = [0, 0, 0, 0]
        negcountsave = np.inf
        for i in range(0, len(idx)-1):
            for ii in [0, idx[0]]:
                for iii in [idx[len(idx)-1], len(x)-1]:
                    for iiii in range(i, len(idx)):
                        z = np.polyfit(
                            x[[ii, idx[i], idx[iiii], iii]], y[[ii, idx[i], idx[iiii], iii]], 3)
                        p = np.poly1d(z)
                        negcount = 0
                        for j in idx:
                            if y[j] - p(x)[j] > 0:
                                negcount += abs(y[j] - p(x)[j])
                            if y[j] - p(x)[j] < 0:
                                if abs(y[j] - p(x)[j]) > 0.01*(max(y)-min(y)):
                                    negcount += np.inf
                                negcount += 40*abs(y[j] - p(x)[j])
                        if z[0] >= 0 and negcount < negcountsave:
                            negcountsave = negcount
                            zsave = z

        baseline = np.poly1d(zsave)(x)
        return baseline

    if(method == 'lin'):
        linesave = np.zeros(len(x))
        anchorpts = np.array([0] + idx.tolist() + [len(x)-1])
        i = 0
        while i != len(anchorpts)-1:
            isave = i
            for j in range(i+1, len(anchorpts)):
                a = (y[anchorpts[j]] - y[anchorpts[i]]) / \
                    (x[anchorpts[j]] - x[anchorpts[i]])
                b = y[anchorpts[i]] - a*x[anchorpts[i]]
                line = a*x[anchorpts[i]:anchorpts[j]+1] + b
                breakcond = False
                for k in range(i, j+1):
                    if y[anchorpts[k]] - line[anchorpts[k]-anchorpts[i]] < -10**-10:
                        breakcond = True
                        break
                if breakcond is False:
                    linesave[anchorpts[i]:anchorpts[j]+1] = line
                    isave = j-1
            i = isave+1
        return linesave


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only
        smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except(ValueError):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(
        -half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def locMin(x, y):
    window = int(0.001*len(x))  # Number of points to average the derivate
    if window <= 1:
        window = 2
    i = 0
    argmin = []
    deriv1, _, _, _, _ = linregress(y[i:i+window], x[i:i+window])
    i += window
    while i < len(x)-window:
        deriv2, _, _, _, _ = linregress(y[i:i+window], x[i:i+window])
        if deriv1 < 0 and deriv2 > 0:
            argmin.append(np.array(y[i-window:i+window]).argmin() + i-window)
        deriv1 = deriv2
        i += window
    return np.array(argmin)
