from dataclasses import dataclass
import numpy as np
from numpy import array as _arr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size']=12

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times"],
    "font.size": 12
})

class Correlation(object):

    def __init__(self, data):

        self.data = _arr(data, dtype=float)
        self.params = {'data':self.data, 'corr_func':None, 'fit_func':None}
        self.results = {}

        self.__set_params()
    
    def __set_params(self):

        def cfunc(v0, vt):
            return np.dot(vt, v0)/np.dot(v0, v0)

        def ffunc(x, *ags):
            return ags[0]*np.power(np.e, x/ags[1])

        self.params.update({'corr_func':cfunc})
        self.params.update({'fit_func':ffunc})

    def analyze(self):

        data = self.params['data']
        cfunc = self.params['corr_func']
        res = []

        for val in data:
            res.append(cfunc(v0=data[0], vt=val))
        self.results.update({'corr_data':res})

        return res

    def fit(self):

        ffunc = self.params['fit_func']

        res = self.results['corr_data']
        xaxis = _arr([i for i in range(len(res))])
        popt, pocv = curve_fit(ffunc, xaxis, res)
        self.results.update({'fit_data':popt})

        return popt

    def plot(self, xlabel, ylabel, xtimes=1, show=True):

        ordata = self.params['data']
        res = self.results['corr_data']
        ffunc = self.params['fit_func']
        popt = self.results['fit_data']

        xaxis = _arr([i for i in range(len(res))])

        fig, ax = plt.subplots(1,1)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
    #     plt.ticklabel_format(style='sci',scilimits=(0,0),axis='both')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        plt.scatter(xaxis, ordata, marker='^', color='#2B547E', facecolor='white')
        plt.plot(xaxis, ffunc(xaxis, *popt), color='red')

        if show:
            plt.show()
        else:
            return plt, ax