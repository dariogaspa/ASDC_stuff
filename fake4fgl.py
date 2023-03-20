import string
import os
#from newsavetxt import *
#import pyfits
from astropy.io import fits #replace pyfits
from numpy import * 
import sys
from math import *
#import numarray
import math
from numpy import *
import numpy as np
import ctx


fgl="/Users/dariogasparrini/Documents/Likelihood_ratio/4LAC-DR4/gll_pscP305uw1410_v1_assoc_v8r7_classes_noex.fits"
num=sys.argv[1]
def sech2(x):
    SECH=2/(exp(x)+exp(-x))
    return pow(SECH,2)

dreq =  3.141592653589793/180.0
dreq1 = 1./(3.1415922653589793/180.0)
hdulist = fits.open(fgl)
tbdata = hdulist[1].data
sourcename=tbdata.field('Source_Name')
name=tbdata.field('NickName')
RA=tbdata.field('RAJ2000')
DEC=tbdata.field('DEJ2000')
LII=tbdata.field('GLON')
BII=tbdata.field('GLAT')
c95sma=tbdata.field('Conf_95_SemiMajor')
c95smi=tbdata.field('Conf_95_SemiMinor')
c95posa=tbdata.field('Conf_95_PosAng')
ts=tbdata.field('Test_Statistic')
pivotene=tbdata.field('Pivot_Energy')
#fluxden=tbdata.field('Flux_Density')
#uncfluxden=tbdata.field('Unc_Flux_Density')
#spi=tbdata.field('Spectral_Index')
#uspi=tbdata.field('Unc_Spectral_Index')
LIIfake=zeros(len(tbdata),float)
BIIfake=zeros(len(tbdata),float)
RAfake=zeros(len(tbdata),float)
DECfake=zeros(len(tbdata),float)
dec = DEC*dreq
ra = RA*dreq
lii=LII*dreq
bii=BII*dreq
b0=5.*dreq
rmax=10.*dreq



for j in range(len(tbdata)):
#for j in range(3):
    galactic=ctx.j20002gal(RA[j],DEC[j])
   # print "RA", "DEC", RA[j], DEC[j],  galactic[0], galactic[1]
    
   # print galactic[0], galactic[1], LII[j], BII[j]
    if abs(BII[j])>=10:
        for k in range(1):
            rad_int=random.randint(1000,5000)
            theta_int=random.randint(0,359999)
            rad=double(rad_int)/1000
            theta=double(theta_int)/1000
            theta_rad=theta*dreq
            RAfake[j]=RA[j]+rad*cos(theta_rad)
            DECfake[j]=DEC[j]+rad*sin(theta_rad)

            if DECfake[j]>89.9:
                print "dai su alto", DECfake[j]
                DECfake[j]=89.9
            if RAfake[j]<0.0:
                print "dai giu de lato sx", RAfake[j]
                RAfake[j]=360.0+RAfake[j]
            if DECfake[j]<-89.9:
                print "dai giu basso", DECfake[j]
                DECfake[j]=-89.9
            if RAfake[j]>360.0:
                print "dai su de lato dx ", RAfake[j]
                RAfake[j]=RAfake[j]-360.0
            galfake=ctx.j20002gal(RAfake[j],DECfake[j])
            #print galfake[0], galfake[1], "             ", j, LII[j], BII[j]              
            LIIfake[j]=galfake[0]
            BIIfake[j]=galfake[1]
               
    else:
         for h in range(1): 
             lintmax=int((lii[j]+(10*dreq))*100000000)
             lintmin=int((lii[j]-(10*dreq))*100000000)
             #print (double(lintmax)/10000)*dreq1, (double(lintmin)/10000)*dreq1, LII[j]
             lii_fake_int=random.randint(lintmin,lintmax)
             #lii_fake_int=random.randint(0,6283)
             lii_fake= double(lii_fake_int)/100000000
             #print lii_fake*dreq1 
             #print lii_fake
             bmax=rmax*(1-sech2(bii[j]/b0))
             if bmax<0.2*dreq:
                 bmax=0.2*dreq
               
                 
             #if lii_fake<lii[j]+(10*dreq) and lii_fake>lii[j]-(10*dreq):
             bintmax=int((bii[j]+bmax)*100000000)
             bintmin=int((bii[j]-bmax)*100000000)
             bmmin=double(bintmin)/100000000
             bmmax=double(bintmax)/100000000
             if bmmax > bmmin:
                 bii_fake_int=random.randint(bintmin,bintmax) 
                 
             else:
                 bii_fake_int=random.randint(bintmax,bintmix)
                #  print "lallo"
             bii_fake= double(bii_fake_int)/100000000   
             if lii_fake<0.0:
                 lii_fake=lii_fake+6.283
             if lii_fake>6.283:
                 lii_fake=lii_fake-6.283
            # print  lii_fake*dreq1, bii_fake*dreq1, j,  LII[j], BII[j], bmax*dreq1 
             #print  lii_fake*dreq1, bii_fake*dreq1,"              ", j, LII[j] , BII[j]
             equatorial=ctx.gal2j2000(lii_fake*dreq1,bii_fake*dreq1)
             LIIfake[j]=lii_fake*dreq1
             BIIfake[j]=bii_fake*dreq1
             RAfake[j]=equatorial[0]
             DECfake[j]=equatorial[1]
             if DECfake[j]>89.9:
                print "dai su alto", DECfake[j]
                DECfake[j]=89.9
	     if DECfake[j]<-89.9:
                print "dai giu basso", DECfake[j]
                DECfake[j]=-89.9
            # if bii_fake<0 and bii_fake>bii[j]-bmax:
            #         print lii_fake*dreq1, bii_fake*dreq1, j
            #         break
            #     if bii_fake>0 and bii_fake<bii[j]+bmax:   
            #         print lii_fake*dreq1, bii_fake*dreq1, j
            #         break
    col0= fits.Column(name='SourceName',format='20A',array=name)
    col1 = fits.Column(name='NickName',format='20A',array=name)
    col2 = fits.Column(name='RA',format='E',array=RAfake)
    col3 = fits.Column(name='DEC',format='E',array=DECfake)
    col4 = fits.Column(name='GLON',format='E',array=LIIfake)
    col5 = fits.Column(name='GLAT',format='E',array=BIIfake)
    col6 = fits.Column(name='Conf_95_SemiMajor',format='E',array=c95sma)
    col7 = fits.Column(name='Conf_95_SemiMinor',format='E',array=c95smi)
    col8 = fits.Column(name='Conf_95_PosAng',format='E',array=c95posa)
    col9 = fits.Column(name='Signif_Avg',format='E',array=ts)
    col10 = fits.Column(name='Pivot_Energy',format='E',array=pivotene)
    #col11 = fits.Column(name='Flux_Density',format='E',array=fluxden)
    #col12 = fits.Column(name='Unc_Flux_Density',format='E',array=uncfluxden)
    #col13 = fits.Column(name='Spectral_Index',format='E',array=spi)
    #col14 = fits.Column(name='Unc_Spectral_Index',format='E',array=uspi)
   
    
    cols=fits.ColDefs([col0,col1, col2, col3, col4, col5,col6, col7, col8, col9, col10]) #col11, col12, col13, col14])
#tbhdu=fits.new_table(cols)
tbhdu=fits.BinTableHDU.from_columns(cols)
nomefits=('LACfake_fake'+num+".fits")
tbhdu.writeto(nomefits)
