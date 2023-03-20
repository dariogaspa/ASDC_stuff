import string
import os
#import pyfits
from astropy.io import fits #replace pyfits
from numpy import * 
import sys
from math import *
#import numarray
import math
from numpy import *
import numpy as np
from matplotlib import pylab
from pylab import load
#from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord #replace PyAstronomy


#hdulistp = pyfits.open('fgl_planck30GHz.fits')

hdulist = fits.open('prova.fits')
#hdulist = pyfits.open('nvss_sel.fits')
#hdulist = fits.open('//Users/dariogasparrini/Documents/Likelihood_ratio/4LAC-DR4/NVSS/nvss_DR4.fits')
tbdata = hdulist[1].data
logLR=zeros(len(tbdata),double)
logLR_b=zeros(len(tbdata),double)
LR=zeros(len(tbdata),float)
LR_b=zeros(len(tbdata),float)
rg_match=zeros(len(tbdata),float)
rg_match_b=zeros(len(tbdata),float)
dist=zeros(len(tbdata),float)
area=zeros(len(tbdata),float)  
area_b=zeros(len(tbdata),float)
angle=zeros(len(tbdata),float)
error_radius_planck=zeros(len(tbdata),float)
error_radius_planck_b=zeros(len(tbdata),float)
rel=zeros(len(tbdata),float)
rel_b=zeros(len(tbdata),float)
den =zeros(len(tbdata),float)
#ra_nvss=tbdata.field('RA')
ra_nvss=tbdata.field('RA_nvss')
#dec_nvss=tbdata.field('DEC')
dec_nvss=tbdata.field('DEC_nvss')
fermi_name=tbdata.field('NickName')
#ra_planck=tbdata.field('RAJ2000')
ra_planck=tbdata.field('RA_fermi')
#dec_planck=tbdata.field('DEJ2000')
dec_planck=tbdata.field('DEC_fermi')
c=SkyCoord(ra_nvss,dec_nvss,unit="deg")
c1=SkyCoord(ra_planck,dec_planck,unit="deg")
flux_nvss=tbdata.field('FLUX')
error_ra_nvss=tbdata.field('RA_ERROR')
error_dec_nvss=tbdata.field('DEC_ERROR')
error_radius_fermi_sma2s=tbdata.field('Conf_95_SemiMajor')/2.4477
error_radius_fermi_smi2s=tbdata.field('Conf_95_SemiMinor')/2.4477
error_radius_fermi_sma=error_radius_fermi_sma2s*0.61
error_radius_fermi_smi=error_radius_fermi_smi2s*0.61
error_radius_nvss=tbdata.field('ERROR_RADIUS')
lognlogs=zeros(700000,float)
lognlos_slope=zeros(700000,float)
expected=zeros(700000,float)
expected_b=zeros(700000,float)
radian=57.29578
num=sys.argv[1]
#num='true'
LRtxt=open("LR_norm"+num+".txt",'w')
LRtxt_b=open("LR_norm_b"+num+".txt",'w')
LRvalue=open("LR_value"+num+".txt",'w')
LRvalue_b=open("LR_value_b"+num+".txt",'w')
print len(tbdata)


for k in range(len(tbdata)):
    error_radius_planck[k]=sqrt(error_radius_fermi_sma[k]*error_radius_fermi_smi[k])
    error_radius_planck_b[k]=sqrt(error_radius_fermi_sma2s[k]*error_radius_fermi_smi2s[k])
    angle[k]=dec_planck[k]/radian
#    dist[k]=sqrt(pow((ra_nvss[k]-ra_planck[k])*cos(angle[k]),2)+pow((dec_nvss[k]-dec_planck[k]),2))*60.0
    dist[k]= c[k].separation(c1[k]).arcmin #usually in hours, converted in arcmin 
   #dist[k]=pyasl.getAngDist(ra_planck[k],dec_planck[k],ra_nvss[k],dec_nvss[k])*60.0
    rg_match[k]=dist[k]/sqrt(pow(error_radius_planck[k]*60.0,2)+pow(error_radius_nvss[k]/60.0,2))
    rg_match_b[k]=dist[k]/sqrt(pow(error_radius_planck_b[k]*60.0,2)+pow(error_radius_nvss[k]/60.0,2))
    den[k]=sqrt(pow(error_radius_planck[k]*60.0,2)+pow(error_radius_nvss[k]/60.0,2))
    if (rg_match[k] > 12.) :
        rg_match[k]=12.
    area[k]=3.14159*error_radius_fermi_sma[k]*error_radius_fermi_smi[k]
    area_b[k]=3.14159*error_radius_fermi_sma2s[k]*error_radius_fermi_smi2s[k]
    if (flux_nvss[k] > 1.0) and (flux_nvss[k] < 6.0):
        lognlogs[k]=0.68918365
        lognlos_slope[k]=-0.70398813
    if (flux_nvss[k] > 6.0) and (flux_nvss[k] < 40.0):
        lognlogs[k]=0.33829486
        lognlos_slope[k]=-0.8440029
    if (flux_nvss[k] > 40.0) and (flux_nvss[k] < 250.0):
        lognlogs[k]=0.10158369
        lognlos_slope[k]=-1.2
    if (flux_nvss[k] > 250.0):
        lognlogs[k]=0.048373964
        lognlos_slope[k]=-1.72
    expected[k]=lognlogs[k]*pow(flux_nvss[k]/1000.0,lognlos_slope[k])*area[k]
    expected_b[k]=lognlogs[k]*pow(flux_nvss[k]/1000.0,lognlos_slope[k])*area_b[k]
    LR[k]=exp(-pow(rg_match[k],2.)/2.)/expected[k]
    LR_b[k]=exp(-pow(rg_match_b[k],2.)/2.)/expected_b[k]
    logLR[k]=log10(LR[k])
    logLR_b[k]=log10(LR_b[k])
#### RISULTATO FIt######
#    p0_b=0.559416
#    p1_b=0.796518
#old rel
    #p0_b=1.1359
    #p1_b=0.8365
#new rel
    p0_b=1.9128871763905912
    p1_b=1.167011122658787
    rel_b[k]=((((1-p0_b*(1/exp(logLR_b[k]*p1_b))))))

    #p0=0.752955
    #p1=0.707942

    p0=1.2640
    p1=0.9869

    rel[k]=((((1-p0*(1/exp(logLR[k]*p1))))))

    print k

#    print ra_nvss[k], dec_nvss[k], logLR[k],logLR_b[k],  area[k],area_b[k], dist[k], rg_match[k],rg_match_b[k], flux_nvss[k],rel[k]


logLR_nz=logLR!=0
ind_logLR=logLR[logLR_nz]
fig1 = pylab.figure(figsize=(10, 10))
LR_hist=pylab.hist(ind_logLR, bins=200, range=(-30,10),  density=False, histtype='step')


for k in range(200):
    LRtxt.write(str(LR_hist[1][k])) 
    LRtxt.write("  ")
    LRtxt.write(str(LR_hist[0][k]))
    LRtxt.write("\n")

for k in range(len(tbdata)):
    LRvalue.write(str(logLR[k])) 
    LRvalue.write("\n")

pylab.savefig("LR"+num+".png")



logLR_b_nz=logLR_b!=0
print "no_INF", logLR_b_nz
ind_logLR_b=logLR_b[logLR_b_nz]
fig2 = pylab.figure(figsize=(10, 10))
LR_hist_b=pylab.hist(ind_logLR_b, bins=200, range=(-30,10),  density=False, histtype='step')


for kk in range(200):
    LRtxt_b.write(str(LR_hist_b[1][kk]))
    LRtxt_b.write("  ")
    LRtxt_b.write(str(LR_hist_b[0][kk]))
    LRtxt_b.write("\n")

for kk in range(len(tbdata)):
    LRvalue_b.write(str(logLR_b[kk]))
    LRvalue_b.write("\n")

pylab.savefig("LR_b"+num+".png")


#####SCOMMENTA DOPO IL FIT DELLA REL#######


#col1 = fits.Column(name='Fermi_name',format='20A',array=fermi_name)
#col2 = fits.Column(name='RA_radio',format='E',array=ra_nvss)
#col3 = fits.Column(name='DEC_radio',format='E',array=dec_nvss)
#col4 = fits.Column(name='Flux',format='E',array=flux_nvss)
#col5 = fits.Column(name='LR',format='E',array=logLR)
#col6 = fits.Column(name='rg_match',format='E',array=rg_match)
#col7 = fits.Column(name='area',format='E',array=area)
#col8 = fits.Column(name='dista',format='E',array=dist)
#col9 = fits.Column(name='LR_b',format='E',array=logLR_b)
#col10 = fits.Column(name='rg_match_b',format='E',array=rg_match_b)
#col11 = fits.Column(name='area_b',format='E',array=area_b)
#col12 = fits.Column(name='rel_b',format='E',array=rel_b)
#col13 = fits.Column(name='rel',format='E',array=rel)
#col14 = fits.Column(name='fermi error ave',format='E',array=error_radius_planck)
#col15 = fits.Column(name='denom',format='E',array=den)
#cols= fits.ColDefs([col1, col2, col3, col4, col5, col6,col7,col8,col9,col10,col11, col12,col13,col14,col15])
#tbhdu=fits.BinTableHDU.from_columns(cols)
#tbhdu.writeto('table_NVSS_DR4.fits')




