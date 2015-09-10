import os 
import sys
import gt_apps as my_apps
import datetime  
#from make1FGLxml import *
import pyLikelihood
from UnbinnedAnalysis import *
from UnbinnedAnalysis import *
from UpperLimits import UpperLimits
import xml_creator_P8_v1
#from likeSED import *
import fileinput
import pyfits
from math import *
from numpy import *


def computeDate(MET):
    
    metdate  = datetime.datetime(2001, 1, 1,0,0,0)
    dt=datetime.timedelta(seconds=MET)
    grb_date=metdate + dt
    yy=grb_date.year
    mm=grb_date.month
    dd=grb_date.day
    hr=grb_date.hour
    mi=grb_date.minute
    ss=grb_date.second
#    fff=float(ss+60.*mi+3600.*hr)/86.4
    date1    = datetime.datetime(yy,mm, dd, hr , mi , ss)
    metdate  = datetime.datetime(2001, 1, 1,0,0,0)
    difference=date1-metdate
    second=float((difference).seconds/86400.)
    MJD= 51910.+float((difference).days)+second
    return MJD



#def sedtool(sedfile):
#  hdulist = pyfits.open(sedfile)
#  tbdata = hdulist[1].data
#  tbdata2 = hdulist[2].data
#  lb=tbdata2.field('Bin Low Edge')
#  hb=tbdata2.field('Bin High Edge')
#  ce=tbdata.field('Center Energy')
#  E2E=tbdata.field('E^2dN/dE')
#  E2E_err=tbdata.field('E^2dN/dE_Error')
#  sedtxt=open("SED.txt",'w')
#  hz=zeros(len(ce),float)
#  freq_err=zeros(len(ce),float)
#  for k in range(len(tbdata)):
#    freq_err[k]=(hb[k]-lb[k])*2.418E+23/2
#    hz[k]=lb[k]*2.418E+23+freq_err[k]
#    if E2E_err[k]!=0:
#       sedtxt.write(str(hz[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(freq_err[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(E2E[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(E2E_err[k]))
#       sedtxt.write("\n")
#    if E2E_err[k]==0:
#       sedtxt.write(str(hz[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(freq_err[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(E2E[k]))
#       sedtxt.write("  ")
#       sedtxt.write(str(E2E_err[k]))
#       sedtxt.write(" UL \n")
     
 


def search(stringsearch,filename):
  for line in fileinput.input(filename):
    if line.find(stringsearch) >= 0:
        lineae=fileinput.filelineno()
        return lineae

kVersion=('\nThis script is simply a wrap-up of commands for Fermi Analyis, version P8 1.0 \n'
          'For info, comments or suggestions: gasparrini@asdc.asi.it and cutini@asdc.asi.it \n')
    
kSynopsis=('\nSYNOPSIS: analysis_fermi.py Source_Name Ra Dec Tmin Tmax EMin Emax SC_file ROI\n'
           'EXAMPLE: analysis_fermi.py 3C279 194.04 -5.789 247104000 248832000 100 100000 L110624055608E0D2F37E81_SC00.fits 10')

def main(NAME,RA,DEC,TSTART,TSTOP,EMIN,EMAX,Np, path, ROIu):
    #outdir = os.environ["FERMI_TMPLATAREA"]
    gtliketxt=open("%s/%s_gtlike.txt"%(path,Np),'w')
    gtsedtxt=open("%s/%s_sed.txt"%(path,Np),'w')
    SCC='%s_SC00.fits'%(Np)
    SC=path+SCC
    Npp=path+Np
    print SC
    ROIue=float(ROIu)+10
    os.system("ls -1 '"+Npp+"'_PH*.fits > %s/%s_events.list" %(path,Np))
 
   # os.system('ls -1 'Np'+'PH*.fits > %s/%s_events.list' %(path,Np)
#GTSELECT
    my_apps.filter['evclass'] = 128
    my_apps.filter['evtype'] = 3
#    my_apps.filter['evclsmin'] = 3
#    my_apps.filter['evclsmax'] = 4
    my_apps.filter['ra'] = RA
    my_apps.filter['dec'] = DEC
    my_apps.filter['rad'] = ROIu
    my_apps.filter['emin'] = EMIN
    my_apps.filter['emax'] = EMAX
    my_apps.filter['zmax'] = 90
    my_apps.filter['tmin'] = TSTART
    my_apps.filter['tmax'] = TSTOP
    my_apps.filter['infile'] = '@%s/%s_events.list' %(path,Np)
    my_apps.filter['outfile'] = '%s/%s_filtered.fits'%(path,Np)
    my_apps.filter.run()
#GTMKTIME
    my_apps.maketime['scfile'] = SC
    my_apps.maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
    my_apps.maketime['roicut'] = 'no'
    my_apps.maketime['evfile'] = '%s/%s_filtered.fits' %(path,Np)
    my_apps.maketime['outfile'] = '%s/%s_filtered_gti.fits' %(path,Np)
    my_apps.maketime.run()
#
#    my_apps.counts_map['evfile'] = '%s/%s_filtered_gti.fits'%(path,Np)
#    my_apps.counts_map['scfile'] = SC
#    my_apps.counts_map['outfile'] = '%s/%s_CountMap.fits'%(path,Np)
#    my_apps.counts_map.run()
#
    my_apps.expCube['evfile'] =  '%s/%s_filtered_gti.fits'%(path,Np)
    my_apps.expCube['scfile'] = SC
    my_apps.expCube['zmax'] = 90
    my_apps.expCube['outfile'] = '%s/%s_expCube.fits' %(path,Np)
    my_apps.expCube['dcostheta'] = 0.025
    my_apps.expCube['binsz'] = 1
    my_apps.expCube.run()

    my_apps.expMap['evfile'] = '%s/%s_filtered_gti.fits'%(path,Np)
    my_apps.expMap['scfile'] = SC
    my_apps.expMap['expcube'] ='%s/%s_expCube.fits'  %(path,Np)
    my_apps.expMap['outfile'] ='%s/%s_expMap.fits'  %(path,Np)
#    my_apps.expMap['irfs'] ='P7REP_SOURCE_V15'
    my_apps.expMap['irfs'] ='CALDB'
    my_apps.expMap['srcrad'] = ROIue
    my_apps.expMap['nlong'] =120
    my_apps.expMap['nlat'] =120
    my_apps.expMap['nenergies'] =20
    my_apps.expMap.run()

    #sara xml model
    roiname='%s/%s_filtered_gti.fits' %(path,Np)
    xml_creator_P8_v1.main(path,NAME,float(RA),float(DEC),float(EMIN), float(EMAX), 20,Np)
    xmlmodelname='%s/%s_model.xml' %(path,Np)
    
    my_apps.diffResps['evfile'] = '%s/%s_filtered_gti.fits'%(path,Np)
    my_apps.diffResps['scfile'] = SC
    my_apps.diffResps['srcmdl'] = xmlmodelname
    my_apps.diffResps['irfs'] = 'CALDB'
    my_apps.diffResps.run()
    
    
    xmlfitname='%s/%s_fit1.xml' %(path,Np)
    expMapFile='%s/%s_expMap.fits'  %(path,Np)
    expCubeFile='%s/%s_expCube.fits'  %(path,Np)
    obs = UnbinnedObs(roiname,SC ,expMap=expMapFile,expCube=expCubeFile,irfs='CALDB')
    like1 = UnbinnedAnalysis(obs,xmlmodelname,optimizer='NewMinuit')
    like1.fit(verbosity=0,optObject=likeobj)
    print likeobj.getRetCode()
    sourceDetails = {}
    for source in like1.sourceNames():
		sourceDetails[source] = like1.Ts(source)
	for source,TS in sourceDetails.iteritems():
		if (TS < 2):
	      print "Deleting...", source, " TS = ", TS
	      like1.deleteSource(source)	
	like1.fit(verbosity=0,optObject=likeobj)
    like1.logLike.writeXml(xmlfitname)

    


 #   numl=search(NAME,xmlfitname)
 #   numlg=str(numl+3)
 #   os.system("sed '"+numlg+","+numlg+" s/free=\"1\"/free=\"0\"/' "+xmlfitname+ " > xml_sed.xml ")
 #   inputs=likeInput(like1,NAME,model="xml_sed.xml",nbins=9,phCorr=1.0)
    #low_edges = [200.,914.61,1955.87,8944.27,19127.05,40902.61]
    #high_edges = [427.69,1955.87,8944.27,19127.05,40902.61,187049.69]
    #centers = [0.2767, 1.265,  5.787, 12.37, 26.46, 86.60]
    #inputs.customBins(low_edges,high_edges)
 #   inputs.plotBins()
 #   inputs.fullFit(CoVar=True)
 #   sed = likeSED(inputs)
 #   sed.getECent()
  #  sed.fitBands()
   # sed.Plot()
    result=like1.model[NAME] 
    TS=like1.Ts(NAME)
#    I = like1.model[NAME].funcs['Spectrum'].getParam('Integral').value()
    flux = like1.flux(NAME,emin=100)  
#    flux=I*1e-9
    
    gamma = like1.model[NAME].funcs['Spectrum'].getParam('Index').value()
    cov_gg =like1.model[NAME].funcs['Spectrum'].getParam('Index').error()
 #   cov_II = like1.model[NAME].funcs['Spectrum'].getParam('Integral').error()
    flux_err = like1.fluxError(NAME,emin=100)
#    flux_err=cov_II*1e-9

    e=1000.0
    a=1
    b=1.e-18
    lenergy_bin=log10(double(EMIN))+(log10(double(EMAX))-log10(double(EMIN)))/2
    energy_bin=pow(10,lenergy_bin)
    freq=2.42e22*energy_bin/100.0
    ums = 1.-gamma
    conv=ums*pow(energy_bin,(-gamma))/(pow(double(EMAX),ums)-pow(double(EMIN),ums))*6.62e-2*(energy_bin/100.0)
 # conv is in Jy
 # now convert in nufnu erg/cm2/s
    convjy=conv*freq*1.e-23
    nufnu=flux*convjy
    b=flux_err*convjy
    err_log=log10((nufnu+b)/nufnu)
     #cout<<freq<<" "<<a<<" "<<nufnu<<" "<<b<<endl;
     #cout<<log10(freq)<<" "<<log10(a)<<" "<<log10(nufnu)<<" "<<err_log<<endl;

    date_start=computeDate(float(TSTART))
    date_stop=computeDate(float(TSTOP))


  #  like1.plot()
  #  fitsedname='%s_9bins_likeSEDout.fits' %NAME
  #  sedtool(fitsedname)


    print NAME, " TS=", TS
#    print result
#    print like1.model
    print "spectral index= ", gamma, " +/-", cov_gg 
    print " Flux=", flux, "+/-", flux_err
    print "freq", freq, " nuFnu=", nufnu, b,
 #   print "'UL': ", results_ul, err
    gtliketxt.write(NAME)
    gtliketxt.write(" RA=")
    gtliketxt.write(RA)
    gtliketxt.write(" DEC= ")
    gtliketxt.write(DEC)
    gtliketxt.write(" TS= ")
    gtliketxt.write(str(TS))
    gtliketxt.write("\n")
    gtliketxt.write(" Time Interval (MJD) ")
    gtliketxt.write(str(date_start))
    gtliketxt.write(" ")
    gtliketxt.write(str(date_stop))
    gtliketxt.write("\n ")
    gtliketxt.write("Flux ")

    if TS <25:
           obs = UnbinnedObs(roiname,SC ,expMap=expMapFile,expCube=expCubeFile,irfs='CALDB')	   
	   like1 = UnbinnedAnalysis(obs,xmlmodelname,optimizer='NewMinuit')
           like1.fit(verbosity=0)
           ul=UpperLimits(like1)
           UL=ul[NAME].compute(emin=double(EMIN),emax=double(EMAX))
           results_ul=UL[1]*1E-9
           err=0
           print "'UL': ", results_ul, err
           gamma_ul=2.0
	   ums_ul = 1.-gamma_ul
           conv_ul=ums_ul*pow(energy_bin,(-gamma_ul))/(pow(double(EMAX),ums_ul)-pow(double(EMIN),ums_ul))*6.62e-2*(energy_bin/100.0)
           convjy_ul=conv_ul*freq*1.e-23
           nufnu_ul=results_ul*convjy_ul
           b=err*convjy_ul
           print "freq", freq,  "0  nuFnu=", nufnu_ul, b,
           gtliketxt.write(str(results_ul)) 
           gtliketxt.write(" 0  ")
           #gtliketxt.write(err)
           gtsedtxt.write(str(freq))
           gtsedtxt.write(" | 0 ")
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(nufnu_ul))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(b))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(date_start))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(date_stop))
           gtsedtxt.write(" | ")
           gtsedtxt.write(" UL ")
           gtsedtxt.write(" | ")

    else:
    
           gtliketxt.write(str(flux))
           gtliketxt.write(" ")
           gtliketxt.write(str(flux_err))
           gtliketxt.write("\n")
           gtliketxt.write("Spectral Index = ")
           gtliketxt.write(str(gamma))
           gtliketxt.write(" ")
           gtliketxt.write(str(cov_gg))
           gtsedtxt.write(" ")
           gtsedtxt.write(str(freq))
           gtsedtxt.write(" | 0 ")
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(nufnu))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(b))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(date_start))
           gtsedtxt.write(" | ")
           gtsedtxt.write(str(date_stop))
           gtsedtxt.write(" | ")


class GenericErr(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)


if __name__ == "__main__":
    argc=len(sys.argv)
    argv=sys.argv
    proceed=True
    #print argc,argv
    if argc==1 or argc<11 or '--help'==argv[1].lower():
        print kVersion+kSynopsis
        proceed=False
    elif '--version'==argv[1].lower():
        print kVersion
        proceed=False
    else:
        NAME=argv[1]
        RA=argv[2]
        DEC=argv[3]
        TSTART=argv[4]
        TSTOP=argv[5]
        EMIN=argv[6]
        EMAX=argv[7]
        Np=argv[8]
        path=argv[9]
        ROIu=argv[10]
        st=os.getenv('FERMI_DIR',0)
     #   print st
        if st==0:
          raise Exception( 'Please initialize Science tools.Can\'t proceed.\n')
        if proceed:
		try: 
			main(NAME,RA,DEC,TSTART,TSTOP,EMIN,EMAX,Np,path,ROIu)
		except RuntimeError as detail:
			errorlogtxt=open("%s/%s_sed_error.txt"%(path,Np),'w')
			errorlogtxt.write(str(detail))
			print detail
        print 'this is the end'
