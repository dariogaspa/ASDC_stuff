import os
import sys
import gt_apps as my_apps
# from make1FGLxml import *
import pyLikelihood
from UnbinnedAnalysis import *
from UpperLimits import UpperLimits
import xml_creator_P8_v1
from likeSED import *
import fileinput
import pyfits
from math import *
from numpy import *


def sedtool(sedfile):
    hdulist = pyfits.open(sedfile)
    tbdata = hdulist[1].data
    tbdata2 = hdulist[2].data
    lb = tbdata2.field('Bin Low Edge')
    hb = tbdata2.field('Bin High Edge')
    ce = tbdata.field('Center Energy')
    E2E = tbdata.field('E^2dN/dE')
    E2E_err = tbdata.field('E^2dN/dE_Error')
    sedtxt = open("SED.txt", 'w')
    hz = zeros(len(ce), float)
    freq_err = zeros(len(ce), float)
    for k in range(len(tbdata)):
        freq_err[k] = (hb[k] - lb[k]) * 2.418E+23 / 2
        hz[k] = lb[k] * 2.418E+23 + freq_err[k]
        if E2E_err[k] != 0:
            sedtxt.write(str(hz[k]))
            sedtxt.write("  ")
            sedtxt.write(str(freq_err[k]))
            sedtxt.write("  ")
            sedtxt.write(str(E2E[k]))
            sedtxt.write("  ")
            sedtxt.write(str(E2E_err[k]))
            sedtxt.write("\n")
        if E2E_err[k] == 0:
            sedtxt.write(str(hz[k]))
            sedtxt.write("  ")
            sedtxt.write(str(freq_err[k]))
            sedtxt.write("  ")
            sedtxt.write(str(E2E[k]))
            sedtxt.write("  ")
            sedtxt.write(str(E2E_err[k]))
            sedtxt.write(" UL \n")


def search(stringsearch, filename):
    for line in fileinput.input(filename):
        if line.find(stringsearch) >= 0:
            lineae = fileinput.filelineno()
            return lineae


kVersion = ('\nThis script is simply a wrap-up of commands for Fermi Analyis, version P8 1.0 \n'
            'For info, comments or suggestions: gasparrini@asdc.asi.it and cutini@asdc.asi.it \n')

kSynopsis = ('\nSYNOPSIS: analysis_fermi_gtlike.py Source_Name Ra Dec Tmin Tmax EMin Emax SC_file ROI\n 1/0'
             'EXAMPLE: analysis_fermi_gtlike.py 3C279 194.04 -5.789 247104000 248832000 100 100000 L110624055608E0D2F37E81_SC00.fits 10 xml=1,0 ')


def main(NAME, RA, DEC, TSTART, TSTOP, EMIN, EMAX, SC, ROIu, xml):
    ROIue = float(ROIu) + 10
    os.system('ls -1 *PH*.fits > %s_events.list' % (NAME))
    my_apps.filter['evclass'] = 128
    my_apps.filter['evtype'] = 3
    my_apps.filter['ra'] = RA
    my_apps.filter['dec'] = DEC
    my_apps.filter['rad'] = ROIu
    my_apps.filter['emin'] = EMIN
    my_apps.filter['emax'] = EMAX
    my_apps.filter['zmax'] = 90
    my_apps.filter['tmin'] = TSTART
    my_apps.filter['tmax'] = TSTOP
    my_apps.filter['infile'] = '@%s_events.list' % (NAME)
    my_apps.filter['outfile'] = '%s_filtered.fits' % (NAME)
    my_apps.filter.run()
    #    maketime
    my_apps.maketime['scfile'] = SC
    my_apps.maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
    my_apps.maketime['roicut'] = 'no'
    my_apps.maketime['evfile'] = '%s_filtered.fits' % (NAME)
    my_apps.maketime['outfile'] = '%s_filtered_gti.fits' % (NAME)
    my_apps.maketime.run()
    #
    my_apps.counts_map['evfile'] = '%s_filtered_gti.fits' % (NAME)
    my_apps.counts_map['scfile'] = SC
    my_apps.counts_map['outfile'] = '%s_CountMap.fits' % (NAME)
    #    my_apps.counts_map.run()
    #
    my_apps.expCube['evfile'] = '%s_filtered_gti.fits' % (NAME)
    my_apps.expCube['scfile'] = SC
    my_apps.expCube['zmax'] = 90
    my_apps.expCube['outfile'] = 'expCube.fits'
    my_apps.expCube['dcostheta'] = 0.025
    my_apps.expCube['binsz'] = 1
    my_apps.expCube.run()

    my_apps.expMap['evfile'] = '%s_filtered_gti.fits' % (NAME)
    my_apps.expMap['scfile'] = SC
    my_apps.expMap['expcube'] = 'expCube.fits'
    my_apps.expMap['outfile'] = 'expMap.fits'
    my_apps.expMap['irfs'] = 'CALDB'
    my_apps.expMap['srcrad'] = ROIue
    my_apps.expMap['nlong'] = 120
    my_apps.expMap['nlat'] = 120
    my_apps.expMap['nenergies'] = 20
    my_apps.expMap.run()

    # sara xml model
    roiname = '%s_filtered_gti.fits' % NAME
    if float(xml) == 0:
        xml_creator_P8_v1.main(NAME, float(RA), float(DEC), float(EMIN), float(EMAX), 15)
        xmlmodelname = '%s_model.xml' % NAME

        my_apps.diffResps['evfile'] = '%s_filtered_gti.fits' % (NAME)
        my_apps.diffResps['scfile'] = SC
        my_apps.diffResps['srcmdl'] = xmlmodelname
        my_apps.diffResps['irfs'] = 'CALDB'
        my_apps.diffResps.run()

        xmlfitname = '%s_fit1.xml' % NAME
        obs = UnbinnedObs(roiname, SC, expMap='expMap.fits', expCube='expCube.fits', irfs='CALDB')
        # like1 = UnbinnedAnalysis(obs,xmlmodelname,optimizer='MINUIT')
        like1 = UnbinnedAnalysis(obs, xmlmodelname, optimizer='NewMinuit')
        likeobj = pyLike.NewMinuit(like1.logLike)
        like1.fit(verbosity=0, optObject=likeobj)
        print likeobj.getRetCode()
        sourceDetails = {}
        for source in like1.sourceNames():
            sourceDetails[source] = like1.Ts(source)
        for source, TS in sourceDetails.iteritems():
            if (TS < 2):
                print "Deleting...", source, " TS = ", TS
                like1.deleteSource(source)
        like1.fit(verbosity=0, optObject=likeobj)
        print "0 is converged", likeobj.getRetCode()
        like1.logLike.writeXml(xmlfitname)

        numl = search(NAME, xmlfitname)
        numlg = str(numl + 3)
        os.system("sed '" + numlg + "," + numlg + " s/free=\"1\"/free=\"0\"/' " + xmlfitname + " > xml_sed.xml ")
        inputs = likeInput(like1, NAME, model="xml_sed.xml", nbins=6, phCorr=1.0)
        inputs.plotBins()
        inputs.fullFit(CoVar=True)
        sed = likeSED(inputs)
        sed.getECent()
        sed.fitBands()
        sed.Plot()
        result = like1.model[NAME]
        TS = like1.Ts(NAME)
        flux = like1.flux(NAME, emin=100)
        gamma = like1.model[NAME].funcs['Spectrum'].getParam('Index').value()
        cov_gg = like1.model[NAME].funcs['Spectrum'].getParam('Index').error()
        #    cov_II = like1.model[NAME].funcs['Spectrum'].getParam('Integral').error()
        flux_err = like1.fluxError(NAME, emin=100)
        like1.plot()
        fitsedname = '%s_6bins_likeSEDout.fits' % NAME
        sedtool(fitsedname)

        print NAME, " TS=", TS
        print result

    if float(xml) == 1:
        xmlmodelname = '%s_model.xml' % NAME
        xmlfitname = '%s_fit1.xml' % NAME
        obs = UnbinnedObs(roiname, SC, expMap='expMap.fits', expCube='expCube.fits', irfs='CALDB')
        # like1 = UnbinnedAnalysis(obs,xmlmodelname,optimizer='MINUIT')
        like1 = UnbinnedAnalysis(obs, xmlmodelname, optimizer='NewMinuit')
        likeobj = pyLike.NewMinuit(like1.logLike)
        like1.fit(verbosity=0, optObject=likeobj)
        print likeobj.getRetCode()
        sourceDetails = {}
        for source in like1.sourceNames():
            sourceDetails[source] = like1.Ts(source)
        for source, TS in sourceDetails.iteritems():
            if (TS < 2):
                print "Deleting...", source, " TS = ", TS
                like1.deleteSource(source)
        like1.fit(verbosity=0, optObject=likeobj)
        print "0 is converged", likeobj.getRetCode()
        like1.logLike.writeXml(xmlfitname)
        numl = search(NAME, xmlfitname)
        numlg = str(numl + 3)
        os.system("sed '" + numlg + "," + numlg + " s/free=\"1\"/free=\"0\"/' " + xmlfitname + " > xml_sed.xml ")
        inputs = likeInput(like1, NAME, model="xml_sed.xml", nbins=6, phCorr=1.0)
        inputs.plotBins()
        inputs.fullFit(CoVar=True)
        sed = likeSED(inputs)
        sed.getECent()
        sed.fitBands()
        sed.Plot()
        result = like1.model[NAME]
        TS = like1.Ts(NAME)
        flux = like1.flux(NAME, emin=100)
        gamma = like1.model[NAME].funcs['Spectrum'].getParam('Index').value()
        cov_gg = like1.model[NAME].funcs['Spectrum'].getParam('Index').error()
        #    cov_II = like1.model[NAME].funcs['Spectrum'].getParam('Integral').error()
        flux_err = like1.fluxError(NAME, emin=100)
        like1.plot()
        fitsedname = '%s_6bins_likeSEDout.fits' % NAME
        sedtool(fitsedname)

        print NAME, " TS=", TS
        print result

        # low_edges = [200.,914.61,1955.87,8944.27,19127.05,40902.61]
        # high_edges = [427.69,1955.87,8944.27,19127.05,40902.61,187049.69]
        # centers = [0.2767, 1.265,  5.787, 12.37, 26.46, 86.60]
        # inputs.customBins(low_edges,high_edges)


class GenericErr(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


if __name__ == "__main__":
    argc = len(sys.argv)
    argv = sys.argv
    proceed = True
    # print argc,argv
    if argc == 1 or argc < 11 or '--help' == argv[1].lower():
        print kVersion + kSynopsis
        proceed = False
    elif '--version' == argv[1].lower():
        print kVersion
        proceed = False
    else:
        NAME = argv[1]
        RA = argv[2]
        DEC = argv[3]
        TSTART = argv[4]
        TSTOP = argv[5]
        EMIN = argv[6]
        EMAX = argv[7]
        SC = argv[8]
        ROIu = argv[9]
        xml = argv[10]
        st = os.getenv('FERMI_DIR', 0)
        #   print st
        if st == 0:
            raise Exception('Please initialize Science tools.Can\'t proceed.\n')
        if proceed: main(NAME, RA, DEC, TSTART, TSTOP, EMIN, EMAX, SC, ROIu, xml)
        print 'this is the end'
