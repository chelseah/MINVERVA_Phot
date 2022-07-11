#!/usr/bin/env python
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import batman 

def get_transit(time, t0, per, rp, a, inc):
    params = batman.TransitParams()
    params.t0 = 0 
    params.per = per
    params.rp = rp
    params.a = a 
    params.inc = inc
    params.ecc = 0.
    params.w = 90.
    params.u = [0.1,0.3]
    params.limb_dark="quadratic"
    m = batman.TransitModel(params, time)
    flux = m.light_curve(params)

    return flux
def bin_lc(x, y, dx=200./60./60./24.):
    minx = np.min(x)
    maxx = np.max(x)
    nbin = int((maxx-minx)/dx)
    binbounds = np.linspace(minx, maxx, nbin)
    binx = np.zeros(len(binbounds)-1)
    biny = np.zeros(len(binbounds)-1)
    binyerr = np.zeros(len(binbounds)-1)
    for i in range(len(binx)):
        index = (x>binbounds[i])*(x<binbounds[i+1])
        binx[i] = np.nanmean(x[index])
        biny[i] = np.nanmean(y[index])
        binyerr[i] = np.nanstd(y[index])/np.sqrt(len(y[index]))

    return [binx, biny, binyerr]

lc1 = pd.read_csv("wasp44_lc.csv", index_col=False)
lc2 = pd.read_csv("wasp36_lc.csv", index_col=False)
lc3 = pd.read_csv("/koi368/xuhuang/MINERVA/MINERVA_phot_data/20180807/wasp44b/wasp44b_measurements.csv", index_col=False)
#for key in lc2.columns:
#    print (key)
#print()
jd1 = np.array(lc1["bjd"])
flux1 = np.array(lc1["12862099"])
flux2 = np.array(lc2["13349647"])
jd2 = np.array(lc2["bjd"])
jd4 = np.array(lc3["BJD_TDB"])
flux4 = np.array(lc3["rel_flux_T1"])/0.0679
binjd2, binflux2, binerr = bin_lc(jd2, flux2)
tc1 = 2458338.102384
p1 = 2.423804
tc2 = 2455569.83771
p2 = 1.53736596
ph1 = (jd1-tc1)%p1
ph1[ph1>p1/2.]-=p1
ph4 = (jd4-tc1)%p1
ph4[ph4>p1/2.]-=p1
ph2 = (jd2-tc2)%p2
ph2[ph2>p2/2.]-=p2
ph3 = (binjd2-tc2)%p2
ph3[ph3>p2/2.]-=p2
#plt.plot(ph4*24, flux4 ,'o', label="2018 CCD WASP-44 r filter transit")
#plt.plot(ph1*24, flux1-0.01, 'o', label="2018 CCD WASP-44 r filter transit")
#plt.plot(ph2/p2, flux2, '.')
#print(binerr)
model = get_transit(ph2, tc2, p2, 0.1367, 5.848, 83.15)

plt.errorbar(ph3*24+0.045, binflux2-0.01, yerr=binerr, fmt='s', label="2021 CMOS WASP-36 transit")
plt.plot(ph2*24, model, '.', color='k')
plt.ylabel("Relative Flux")
plt.xlabel("T-Tc [Hours]")
plt.legend()
plt.show()
