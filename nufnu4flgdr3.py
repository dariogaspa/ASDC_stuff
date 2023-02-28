#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:20:40 2023

@author: dariogasparrini
"""
from astropy.io import fits #replace pyfits
from numpy import zeros,double 
#from matplotlib import pylab

#import matplotlib.pyplot as plt
from math import log10
#from PyAstronomy import pyasl
#from astropy.coordinates import SkyCoord #replace PyAstronomy

hdulist = fits.open('/Users/dariogasparrini/Downloads/gll_psc_v31_exp.fits')
tbdata = hdulist[1].data
orig_cols = tbdata.columns
hdulist.close()

NuFNu_SED_1=zeros(len(tbdata),double)
NuFNu_SED_2=zeros(len(tbdata),double)
NuFNu_SED_3=zeros(len(tbdata),double)
NuFNu_SED_4=zeros(len(tbdata),double)
NuFNu_SED_5=zeros(len(tbdata),double)
NuFNu_SED_6=zeros(len(tbdata),double)
NuFNu_SED_7=zeros(len(tbdata),double)
NuFNu_SED_8=zeros(len(tbdata),double)
NuFNuErrPos_SED_1=zeros(len(tbdata),double)
NuFNuErrPos_SED_2=zeros(len(tbdata),double)
NuFNuErrPos_SED_3=zeros(len(tbdata),double)
NuFNuErrPos_SED_4=zeros(len(tbdata),double)
NuFNuErrPos_SED_5=zeros(len(tbdata),double)
NuFNuErrPos_SED_6=zeros(len(tbdata),double)
NuFNuErrPos_SED_7=zeros(len(tbdata),double)
NuFNuErrPos_SED_8=zeros(len(tbdata),double)
NuFNuErrNeg_SED_1=zeros(len(tbdata),double)
NuFNuErrNeg_SED_2=zeros(len(tbdata),double)
NuFNuErrNeg_SED_3=zeros(len(tbdata),double)
NuFNuErrNeg_SED_4=zeros(len(tbdata),double)
NuFNuErrNeg_SED_5=zeros(len(tbdata),double)
NuFNuErrNeg_SED_6=zeros(len(tbdata),double)
NuFNuErrNeg_SED_7=zeros(len(tbdata),double)
NuFNuErrNeg_SED_8=zeros(len(tbdata),double)

Name=tbdata.field('Source_name')
SqrTs_1=tbdata.field('Sqrt_TS_Band_1')
flux_1=tbdata.field('Flux_Band_1')
errflux_neg_1=tbdata.field('Unc_Flux_Band_1')
errflux_pos_1=tbdata.field('Unc_Flux_Band_2')
nuFnu_1=tbdata.field('nuFnu_Band_1')

SqrTs_2=tbdata.field('Sqrt_TS_Band_2')
flux_2=tbdata.field('Flux_Band_2')
errflux_neg_2=tbdata.field('Unc_Flux_Band_3')
errflux_pos_2=tbdata.field('Unc_Flux_Band_4')
nuFnu_2=tbdata.field('nuFnu_Band_2')

SqrTs_3=tbdata.field('Sqrt_TS_Band_3')
flux_3=tbdata.field('Flux_Band_3')
errflux_neg_3=tbdata.field('Unc_Flux_Band_5')
errflux_pos_3=tbdata.field('Unc_Flux_Band_6')
nuFnu_3=tbdata.field('nuFnu_Band_3')

SqrTs_4=tbdata.field('Sqrt_TS_Band_4')
flux_4=tbdata.field('Flux_Band_4')
errflux_neg_4=tbdata.field('Unc_Flux_Band_7')
errflux_pos_4=tbdata.field('Unc_Flux_Band_8')
nuFnu_4=tbdata.field('nuFnu_Band_4')

SqrTs_5=tbdata.field('Sqrt_TS_Band_5')
flux_5=tbdata.field('Flux_Band_5')
errflux_neg_5=tbdata.field('Unc_Flux_Band_9')
errflux_pos_5=tbdata.field('Unc_Flux_Band_10')
nuFnu_5=tbdata.field('nuFnu_Band_5')

SqrTs_6=tbdata.field('Sqrt_TS_Band_6')
flux_6=tbdata.field('Flux_Band_6')
errflux_neg_6=tbdata.field('Unc_Flux_Band_11')
errflux_pos_6=tbdata.field('Unc_Flux_Band_12')
nuFnu_6=tbdata.field('nuFnu_Band_6')

SqrTs_7=tbdata.field('Sqrt_TS_Band_7')
flux_7=tbdata.field('Flux_Band_7')
errflux_neg_7=tbdata.field('Unc_Flux_Band_13')
errflux_pos_7=tbdata.field('Unc_Flux_Band_14')
nuFnu_7=tbdata.field('nuFnu_Band_7')

SqrTs_8=tbdata.field('Sqrt_TS_Band_8')
flux_8=tbdata.field('Flux_Band_8')
errflux_neg_8=tbdata.field('Unc_Flux_Band_15')
errflux_pos_8=tbdata.field('Unc_Flux_Band_16')
nuFnu_8=tbdata.field('nuFnu_Band_8')

for k in range(len(tbdata)):
     if (SqrTs_1[k] >= 2.) :
         NuFNuErrPos_SED_1[k] = (errflux_pos_1[k]/flux_1[k])*nuFnu_1[k]
         NuFNuErrNeg_SED_1[k] = (errflux_neg_1[k]/flux_1[k])*nuFnu_1[k]
         NuFNu_SED_1[k] = nuFnu_1[k]
     else:
        NuFNuErrPos_SED_1[k] = 0.0
        NuFNuErrNeg_SED_1[k] = 0.0
        NuFNu_SED_1[k]= nuFnu_1[k]*(1+2*(errflux_pos_1[k]/flux_1[k]))
     if (SqrTs_2[k] >= 2.) :
         NuFNuErrPos_SED_2[k] = (errflux_pos_2[k]/flux_2[k])*nuFnu_2[k]
         NuFNuErrNeg_SED_2[k] = (errflux_neg_2[k]/flux_2[k])*nuFnu_2[k]
         NuFNu_SED_2[k] = nuFnu_2[k]
     else:
        NuFNuErrPos_SED_2[k] = 0.0
        NuFNuErrNeg_SED_2[k] = 0.0
        NuFNu_SED_2[k]= nuFnu_2[k]*(1+2*(errflux_pos_2[k]/flux_2[k]))
     if (SqrTs_3[k] >= 2.) :
         NuFNuErrPos_SED_3[k] = (errflux_pos_3[k]/flux_3[k])*nuFnu_3[k]
         NuFNuErrNeg_SED_3[k] = (errflux_neg_3[k]/flux_3[k])*nuFnu_3[k]
         NuFNu_SED_3[k] = nuFnu_3[k]
     else:
        NuFNuErrPos_SED_3[k] = 0.0
        NuFNuErrNeg_SED_3[k] = 0.0
        NuFNu_SED_3[k]= nuFnu_3[k]*(1+2*(errflux_pos_3[k]/flux_3[k]))
     if (SqrTs_4[k] >= 2.) :
         NuFNuErrPos_SED_4[k] = (errflux_pos_4[k]/flux_4[k])*nuFnu_4[k]
         NuFNuErrNeg_SED_4[k] = (errflux_neg_4[k]/flux_4[k])*nuFnu_4[k]
         NuFNu_SED_4[k] = nuFnu_4[k]
     else:
        NuFNuErrPos_SED_4[k] = 0.0
        NuFNuErrNeg_SED_4[k] = 0.0
        NuFNu_SED_4[k]= nuFnu_4[k]*(1+2*(errflux_pos_4[k]/flux_4[k]))
     if (SqrTs_5[k] >= 2.) :
         NuFNuErrPos_SED_5[k] = (errflux_pos_5[k]/flux_5[k])*nuFnu_5[k]
         NuFNuErrNeg_SED_5[k] = (errflux_neg_5[k]/flux_5[k])*nuFnu_5[k]
         NuFNu_SED_5[k] = nuFnu_5[k]
     else:
        NuFNuErrPos_SED_5[k] = 0.0
        NuFNuErrNeg_SED_5[k] = 0.0
        NuFNu_SED_5[k]= nuFnu_5[k]*(1+2*(errflux_pos_5[k]/flux_5[k]))
     if (SqrTs_6[k] >= 2.) :
         NuFNuErrPos_SED_6[k] = (errflux_pos_6[k]/flux_6[k])*nuFnu_6[k]
         NuFNuErrNeg_SED_6[k] = (errflux_neg_6[k]/flux_6[k])*nuFnu_6[k]
         NuFNu_SED_6[k] = nuFnu_6[k]
     else:
        NuFNuErrPos_SED_6[k] = 0.0
        NuFNuErrNeg_SED_6[k] = 0.0
        NuFNu_SED_6[k]= nuFnu_6[k]*(1+2*(errflux_pos_6[k]/flux_6[k]))
     if (SqrTs_7[k] >= 2.) :
         NuFNuErrPos_SED_7[k] = (errflux_pos_7[k]/flux_7[k])*nuFnu_7[k]
         NuFNuErrNeg_SED_7[k] = (errflux_neg_7[k]/flux_7[k])*nuFnu_7[k]
         NuFNu_SED_7[k] = nuFnu_7[k]
     else:
        NuFNuErrPos_SED_7[k] = 0.0
        NuFNuErrNeg_SED_7[k] = 0.0
        NuFNu_SED_7[k]= nuFnu_7[k]*(1+2*(errflux_pos_7[k]/flux_7[k]))
     if (SqrTs_8[k] >= 2.) :
         NuFNuErrPos_SED_8[k] = (errflux_pos_8[k]/flux_8[k])*nuFnu_8[k]
         NuFNuErrNeg_SED_8[k] = (errflux_neg_8[k]/flux_8[k])*nuFnu_8[k]
         NuFNu_SED_8[k] = nuFnu_8[k]
     else:
        NuFNuErrPos_SED_8[k] = 0.0
        NuFNuErrNeg_SED_8[k] = 0.0
        NuFNu_SED_8[k]= nuFnu_8[k]*(1+2*(errflux_pos_8[k]/flux_8[k]))
        
#conversione in GeV  
EGeV=1./1000.        
band_1Min=50.
band_1Max=100.
center_band_1 =10.**((log10(band_1Min)+log10(band_1Max))/2)*(EGeV)
band_2Min=100.
band_2Max=300.
center_band_2 =10.**((log10(band_2Min)+log10(band_2Max))/2)*(EGeV)
band_3Min=300.
band_3Max=1000.
center_band_3 =10.**((log10(band_3Min)+log10(band_3Max))/2)*(EGeV)
band_4Min=1.
band_4Max=3.
center_band_4 =10.**((log10(band_4Min)+log10(band_4Max))/2)
band_5Min=3.
band_5Max=10.
center_band_5 =10.**((log10(band_5Min)+log10(band_5Max))/2)
band_6Min=10.
band_6Max=30.
center_band_6 =10.**((log10(band_6Min)+log10(band_6Max))/2)
band_7Min=30.
band_7Max=100.
center_band_7 =10.**((log10(band_7Min)+log10(band_7Max))/2)
band_8Min=100.
band_8Max=1000.
center_band_8 =10.**((log10(band_8Min)+log10(band_8Max))/2)

#gg=6
#x=[center_band_1, center_band_2, center_band_3,center_band_4,center_band_5, center_band_6, center_band_7, center_band_8  ]
#y=[NuFNu_SED_1[gg],NuFNu_SED_2[gg],NuFNu_SED_3[gg],NuFNu_SED_4[gg],NuFNu_SED_5[gg],NuFNu_SED_6[gg],NuFNu_SED_7[gg],NuFNu_SED_8[gg]]
#y_err=[NuFNuErrPos_SED_1[gg],NuFNuErrPos_SED_2[gg],NuFNuErrPos_SED_3[gg],NuFNuErrPos_SED_4[gg],NuFNuErrPos_SED_5[gg],NuFNuErrPos_SED_6[gg],NuFNuErrPos_SED_7[gg],NuFNuErrPos_SED_8[gg]]
#plt.errorbar(x, y, yerr=y_err, color='g')
#plt.title(Name[gg], fontsize=15)
#plt.loglog(nonposx='clip', nonposy='clip')

#cols = hdulist[1].columns        
       
        
#NuFNu_1 = fits.Column(name='NuFNu_SED_1', format='E',array=NuFNu_SED_1)

#table.add_column(NuFNu_1)

NuFNu_1 = fits.Column(name='NuFNu_SED_1', format='E',array=NuFNu_SED_1)
NuFNu_1_pos = fits.Column(name='NuFNu_SED_1_ErrPos', format='E',array=NuFNuErrPos_SED_1)
NuFNu_1_neg = fits.Column(name='NuFNu_SED_1_ErrNeg', format='E',array=NuFNuErrNeg_SED_1)
NuFNu_2 = fits.Column(name='NuFNu_SED_2', format='E',array=NuFNu_SED_2) 
NuFNu_2_pos = fits.Column(name='NuFNu_SED_2_ErrPos', format='E',array=NuFNuErrPos_SED_2)
NuFNu_2_neg = fits.Column(name='NuFNu_SED_2_ErrNeg', format='E',array=NuFNuErrNeg_SED_2)
NuFNu_3 = fits.Column(name='NuFNu_SED_3', format='E',array=NuFNu_SED_3) 
NuFNu_3_pos = fits.Column(name='NuFNu_SED_3_ErrPos', format='E',array=NuFNuErrPos_SED_3)
NuFNu_3_neg = fits.Column(name='NuFNu_SED_3_ErrNeg', format='E',array=NuFNuErrNeg_SED_3)
NuFNu_4 = fits.Column(name='NuFNu_SED_4', format='E',array=NuFNu_SED_4) 
NuFNu_4_pos = fits.Column(name='NuFNu_SED_4_ErrPos', format='E',array=NuFNuErrPos_SED_4)
NuFNu_4_neg = fits.Column(name='NuFNu_SED_4_ErrNeg', format='E',array=NuFNuErrNeg_SED_4)
NuFNu_5 = fits.Column(name='NuFNu_SED_5', format='E',array=NuFNu_SED_5)
NuFNu_5_pos = fits.Column(name='NuFNu_SED_5_ErrPos', format='E',array=NuFNuErrPos_SED_5)
NuFNu_5_neg = fits.Column(name='NuFNu_SED_5_ErrNeg', format='E',array=NuFNuErrNeg_SED_5)
NuFNu_6 = fits.Column(name='NuFNu_SED_6', format='E',array=NuFNu_SED_6) 
NuFNu_6_pos = fits.Column(name='NuFNu_SED_6_ErrPos', format='E',array=NuFNuErrPos_SED_6)
NuFNu_6_neg = fits.Column(name='NuFNu_SED_6_ErrNeg', format='E',array=NuFNuErrNeg_SED_6)
NuFNu_7 = fits.Column(name='NuFNu_SED_7', format='E',array=NuFNu_SED_7) 
NuFNu_7_pos = fits.Column(name='NuFNu_SED_7_ErrPos', format='E',array=NuFNuErrPos_SED_7)
NuFNu_7_neg = fits.Column(name='NuFNu_SED_7_ErrNeg', format='E',array=NuFNuErrNeg_SED_7)
NuFNu_8 = fits.Column(name='NuFNu_SED_8', format='E',array=NuFNu_SED_8) 
NuFNu_8_pos = fits.Column(name='NuFNu_SED_8_ErrPos', format='E',array=NuFNuErrPos_SED_8)
NuFNu_8_neg = fits.Column(name='NuFNu_SED_8_ErrNeg', format='E',array=NuFNuErrNeg_SED_8)



cols = fits.ColDefs([NuFNu_1, NuFNu_1_pos, NuFNu_1_neg, NuFNu_2, NuFNu_2_pos, NuFNu_2_neg,NuFNu_3, NuFNu_3_pos, NuFNu_3_neg, NuFNu_4, NuFNu_4_pos, NuFNu_4_neg,NuFNu_5, NuFNu_5_pos, NuFNu_5_neg, NuFNu_6, NuFNu_6_pos, NuFNu_6_neg,NuFNu_7, NuFNu_7_pos, NuFNu_7_neg, NuFNu_8, NuFNu_8_pos, NuFNu_8_neg])
tbhdu=fits.BinTableHDU.from_columns(cols)
tbhdu.writeto('only_nufnu_sed.fits')
hdu = fits.BinTableHDU.from_columns(orig_cols + cols)
hdu.writeto('gll_psc_v31_exp_SED.fits')