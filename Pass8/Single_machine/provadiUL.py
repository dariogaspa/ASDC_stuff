import pyLikelihood
from UnbinnedAnalysis import *
from UpperLimits import UpperLimits

def main(NAME,RA,DEC,EMIN,EMAX,SC, ROIfile, xmlmodelname, expmap, expcube):

    obs = UnbinnedObs(ROIfile,SC ,expMap=expmap,expCube=expcube,irfs='P7REP_SOURCE_V15')
    like1 = UnbinnedAnalysis(obs,xmlmodelname,optimizer='MINUIT')
    like1.fit(verbosity=0)
    ul=UpperLimits(like1)
    UL=ul[NAME].compute(emin=float(EMIN),emax=float(EMAX))
    results_ul=UL[1]*1E-9
    err=0
    print "'UL': ", results_ul



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
        EMIN=argv[4]
        EMAX=argv[5]
        SC=argv[6]
        ROIfile=argv[7]
        xmlmodelname=argv[8]
        expmap=argv[9]
        expcube=argv[10]
    # st=os.getenv('FERMI_DIR',0)
     #   print st
      #  if st==0:
      #    raise Exception( 'Please initialize Science tools.Can\'t proceed.\n')
        if proceed: main(NAME,RA,DEC,EMIN,EMAX,SC,ROIfile,xmlmodelname,expmap,expcube)
        print 'this is the end'
