
import string
import os
import pyfits
from numpy import * 
import sys
from math import *
# import numarray
import math
from numpy import *
import numpy as np

def main(path,NAME,RA,DEC,EMIN,EMAX,ROI,Np):
    #bkghome="/usr/local/soft/web/apache-tomcat-6.0.35/apps/fermi/bkg"
    bkghome="/usr/local/soft/web/apache-tomcat-7.0.34/apps/fermi/bkg"
    #bkghome="/home/fermi/prove/bkg/"
    fgl = "%s/gll_psc_v16.fit" %(bkghome)
    dreq =  3.141592653589793/180.0
    dreq1 = 1./(3.1415922653589793/180.0)
    hdulist = pyfits.open(fgl)
    tbdata = hdulist[1].data
    name1fgl = tbdata.field('Source_Name')
    RA1fgl = tbdata.field('RAJ2000')
    DEC1fgl = tbdata.field('DEJ2000')
    LII1fgl = tbdata.field('GLON')
    BII1fgl = tbdata.field('GLAT')
    c68sma = tbdata.field('Conf_68_SemiMajor')
    c68smi = tbdata.field('Conf_68_SemiMinor')
    c68posa = tbdata.field('Conf_68_PosAng')
    c95sma = tbdata.field('Conf_95_SemiMajor')
    c95smi = tbdata.field('Conf_95_SemiMinor')
    c95posa = tbdata.field('Conf_95_PosAng')
    pivotene = tbdata.field('Pivot_Energy')
    fluxden = tbdata.field('Flux_Density')
    uncfluxden = tbdata.field('Unc_Flux_Density')
    spi=tbdata.field('Spectral_Index')
    uspi=tbdata.field('Unc_Spectral_Index')
    flux1000=tbdata.field('Flux1000')
    unflux1000=tbdata.field('Unc_Flux1000')
    vindex=tbdata.field('Variability_Index')
    sigif=tbdata.field('Signif_Peak')
    fluxpeak=tbdata.field('Flux_Peak')
    Vx1fgl=zeros(len(RA1fgl),float)
    Vy1fgl=zeros(len(RA1fgl),float)
    Vz1fgl=zeros(len(RA1fgl),float)
    argx=zeros(len(RA1fgl),float)
    argy=zeros(len(RA1fgl),float)
    argz=zeros(len(RA1fgl),float)
    arg=zeros(len(RA1fgl),float)
    DELTA=zeros(len(RA1fgl),float)
    ra1fgl = RA1fgl*dreq
    dec1fgl = DEC1fgl*dreq
    ra=RA*dreq
    dec=DEC*dreq
    min=10000
    gal="%s/gll_iem_v06.fits" %(bkghome)
    iso="%s/iso_P8R2_SOURCE_V6_v06.txt" %(bkghome)
    nomefile='%s/%s_model.xml' %(path,Np)
    model=open(nomefile,'w')
    model.write('<?xml version="1.0" ?>\n')
    model.write('<source_library title="source library">\n')
    model.write('  <source name="')
    model.write(NAME)
    model.write('" type="PointSource">\n')
    model.write('    <spectrum type="PowerLaw2">\n')
    model.write('      <parameter error="0.00" free="1" max="100000" min="1e-07" name="Integral" scale="1e-09" value="1.000000"/>\n')
    model.write('      <parameter error="0.00" free="1" max="5" min="-1" name="Index" scale="-1" value="2.000"/>\n')
    model.write('     <parameter free="0" max="3e6" min="20" name="LowerLimit" scale="1" value="')
    model.write(str(EMIN))
    model.write('"/>\n')
    model.write('      <parameter free="0" max="3e6" min="20" name="UpperLimit" scale="1" value="')
    model.write(str(EMAX))
    model.write('"/>\n')
    model.write('    </spectrum>\n')
    model.write('    <spatialModel type="SkyDirFunction">\n')
    model.write('      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="')
    model.write(str(RA))
    model.write('"/>\n')
    model.write('      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="')
    model.write(str(DEC))
    model.write('"/>\n')
    model.write('    </spatialModel>\n')
    model.write('  </source>\n')
    for k in range(len(tbdata)):
   #     Vx = cos(dec)*cos(ra)
   #     Vy = cos(dec)*sin(ra)
   #     Vz = cos(dec)
   #     Vx1fgl[k] = cos(dec1fgl[k])*cos(ra1fgl[k])
   #     Vy1fgl[k] = cos(dec1fgl[k])*sin(ra1fgl[k])
   #     Vz1fgl[k] = sin(dec1fgl[k])
   #     argx[k] = Vx1fgl[k]*Vx
   #     argy[k] = Vy1fgl[k]*Vy
   #     argz[k] = Vz1fgl[k]*Vz
   #     arg[k] = (argx[k]+argy[k]+argz[k])#-0.0000001
   #     if arg[k]>=1:
   #         DELTA[k] = 0.000000000001
    #    if  arg[k]<1:
    #        DELTA[k] = acos(arg[k])*dreq1
        DELTA[k]=sqrt(pow((RA-RA1fgl[k])*cos(dec),2)+pow((DEC-DEC1fgl[k]),2))
        if  DELTA[k]<min:
            min=DELTA[k]
            namemin=name1fgl[k]
            c95smamin=c95sma[k]
            #print  DELTA[k], name1fgl[k]
            print min/c95smamin, namemin, c95smamin, RA1fgl[k], DEC1fgl[k],  DELTA[k] 
        if  min/c95smamin<1:
            if DELTA[k] < 10. and DELTA[k]>min :
                model.write('<!-- distance =')
                model.write(str(DELTA[k]))
                model.write('-->\n')
                model.write('  <source name="')
                model.write(name1fgl[k])
                model.write('" type="PointSource">\n')
                model.write('    <spectrum type="PowerLaw2">\n')
                model.write('      <parameter error="0.00" free="1" max="1000" min="1e-09" name="Integral" scale="1e-06" value="')
                model.write(str(flux1000[k]/1e-06))
                model.write('"/>\n')
                model.write('      <parameter error="0.00" free="0" max="5" min="1" name="Index" scale="-1" value="')
                model.write(str(spi[k]))
                model.write('"/>\n')
                model.write('     <parameter free="0" max="3e6" min="20" name="LowerLimit" scale="1" value="')
                model.write(str(EMIN))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="3e6" min="20" name="UpperLimit" scale="1" value="')
                model.write(str(EMAX))
                model.write('"/>\n')
                model.write('    </spectrum>\n')
                model.write('    <spatialModel type="SkyDirFunction">\n')
                model.write('      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="')
                model.write(str(RA1fgl[k]))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="')
                model.write(str(DEC1fgl[k]))
                model.write('"/>\n')
                model.write('    </spatialModel>\n')
                model.write('  </source>\n')
                
            if DELTA[k] > 10. and DELTA[k] < float(ROI):
                    model.write('<!-- distance =')
                    model.write(str(DELTA[k]))
                    model.write('-->\n')    
                    model.write('  <source name="')
                    model.write(name1fgl[k])
                    model.write('" type="PointSource">\n')
                    model.write('    <spectrum type="PowerLaw2">\n')
                    model.write('      <parameter error="0.00" free="0" max="1000" min="1e-09" name="Integral" scale="1e-06" value="')
                    model.write(str(flux1000[k]/1e-06))
                    model.write('"/>\n')
                    model.write('      <parameter error="0.00" free="0" max="5" min="0" name="Index" scale="-1" value="')
                    model.write(str(spi[k]))
                    model.write('"/>\n')
                    model.write('     <parameter free="0" max="3e6" min="20" name="LowerLimit" scale="1" value="')
                    model.write(str(EMIN))
                    model.write('"/>\n')
                    model.write('      <parameter free="0" max="3e6" min="20" name="UpperLimit" scale="1" value="')
                    model.write(str(EMAX))
                    model.write('"/>\n')
                    model.write('    </spectrum>\n')
                    model.write('    <spatialModel type="SkyDirFunction">\n')
                    model.write('      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="')
                    model.write(str(RA1fgl[k]))
                    model.write('"/>\n')
                    model.write('      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="')
                    model.write(str(DEC1fgl[k]))
                    model.write('"/>\n')
                    model.write('    </spatialModel>\n')
                    model.write('  </source>\n')
        if  min/c95smamin>1:
            if DELTA[k] < 10.:
                model.write('<!-- distance =')
                model.write(str(DELTA[k]))
                model.write('-->\n')
                model.write('  <source name="')
                model.write(name1fgl[k])
                model.write('" type="PointSource">\n')
                model.write('    <spectrum type="PowerLaw2">\n')
                model.write('      <parameter error="0.00" free="1" max="1000" min="1e-09" name="Integral" scale="1e-06" value="')
                model.write(str(flux1000[k]/1e-6))
                model.write('"/>\n')
                model.write('      <parameter error="0.00" free="0" max="5" min="0" name="Index" scale="-1" value="')
                model.write(str(spi[k]))
                model.write('"/>\n')
                model.write('     <parameter free="0" max="3e6" min="20" name="LowerLimit" scale="1" value="')
                model.write(str(EMIN))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="3e6" min="20" name="UpperLimit" scale="1" value="')
                model.write(str(EMAX))
                model.write('"/>\n')
                model.write('    </spectrum>\n')
                model.write('    <spatialModel type="SkyDirFunction">\n')
                model.write('      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="')
                model.write(str(RA1fgl[k]))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="')
                model.write(str(DEC1fgl[k]))
                model.write('"/>\n')
                model.write('    </spatialModel>\n')
                model.write('  </source>\n')            
            if DELTA[k] > 10. and DELTA[k] < 20:
                model.write('<!-- distance =')
                model.write(str(DELTA[k]))
                model.write('-->\n')    
                model.write('  <source name="')
                model.write(name1fgl[k])
                model.write('" type="PointSource">\n')
                model.write('    <spectrum type="PowerLaw2">\n')
                model.write('      <parameter error="0.00" free="0" max="1000" min="1e-09" name="Integral" scale="1e-06" value="')
                model.write(str(flux1000[k]/1e-06))
                model.write('"/>\n')
                model.write('      <parameter error="0.00" free="0" max="5" min="0" name="Index" scale="-1" value="')
                model.write(str(spi[k]))
                model.write('"/>\n')
                model.write('     <parameter free="0" max="3e6" min="20" name="LowerLimit" scale="1" value="')
                model.write(str(EMIN))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="3e6" min="20" name="UpperLimit" scale="1" value="')
                model.write(str(EMAX))
                model.write('"/>\n')
                model.write('    </spectrum>\n')
                model.write('    <spatialModel type="SkyDirFunction">\n')
                model.write('      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="')
                model.write(str(RA1fgl[k]))
                model.write('"/>\n')
                model.write('      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="')
                model.write(str(DEC1fgl[k]))
                model.write('"/>\n')
                model.write('    </spatialModel>\n')
                model.write('  </source>\n')
# insert the Diffuse Model
  #  model.write('<source name="gll_iem_v05" type="DiffuseSource">\n')
  #  model.write('   <spectrum type="PowerLaw">\n')
  #  model.write('	<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1."/>\n')
  #  model.write('      <parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0."/>\n')
  #  model.write('      <parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n')
  #  model.write('   </spectrum>\n')
  #  model.write('   <spatialModel file="')
  #  model.write(str(gal))
  #  model.write('" type="MapCubeFunction">\n')
  #  model.write('      <parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n')
  #  model.write('   </spatialModel>\n')
  #  model.write('  </source>\n')

    model.write('<source name="gll_iem_v06" type="DiffuseSource">\n')
    model.write('<spectrum type="PowerLaw">\n')
    model.write('<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1"/>\n')
    model.write('<parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n')
    model.write('<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n')
    model.write('</spectrum>\n')
    model.write('   <spatialModel file="')
    model.write(str(gal))
    model.write('" type="MapCubeFunction">\n')
    model.write('<parameter free="0" max="1000.0" min="0.001" name="Normalization" scale= "1.0" value="1.0"/>\n')
    model.write('</spatialModel>\n')
    model.write('</source>\n')

    model.write('<source name="iso_P8R2_SOURCE_V6_v06" type="DiffuseSource">\n')
    model.write('   <spectrum file="')
    model.write(str(iso))
    model.write('" type="FileFunction">\n')
    model.write('    <parameter free="1" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n')
    model.write('   </spectrum>\n')
    model.write('   <spatialModel type="ConstantValue">\n')
    model.write('    <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n')
    model.write('   </spatialModel>\n')
    model.write('  </source>\n')
    model.write('</source_library>\n')
      
    
if __name__ == "__main__":
    argc = len(sys.argv)    
    if argc !=8:
        print "python xml editor  using 1FGL catalog"
        print "insert path NAME_SOURCE RA DEC EMIN EMAX ROI Np"
        print "troubles contact: sara.cutini@asdc.asi.it"
    if argc==7:
        path=sys.argv[1]
        NAME=sys.argv[2]
        RA=float(sys.argv[3])
        DEC=float(sys.argv[4])
        EMIN=sys.argv[5]
        EMAX=sys.argv[6]
        ROI=sys.argv[7]
        Np=sys.argv[8]
        main(path,NAME,RA,DEC,EMIN,EMAX,ROI,Np)
