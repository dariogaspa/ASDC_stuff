#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
#include<stdio.h>
#include<stdlib.h>
 
using namespace std;

int main(int argc, char *argv[]) {
 double freq, emin, emax,ee,aaa;
  FILE *fp;
  double dato;
  double time[40000];
  double UL[40000];
  double tstop[40000];
  double tstart[40000];
  double delta[40000];
  double flux[40000];
  double fluxerr[40000];
  double ra[40000];
  double dec[40000];
  double a[40000][4];
  double flag[40000];
  double sp[40000];
  double mjdsta[40000];
  double mjdsto[40000];
  double sperr[40000];
  double b[40000];
  double convjy[40000];
  double conv[40000];
  double ums[40000];
  double nufnu[40000];
  int nbins=395;
  // cin>>nbins;

  fp = fopen(argv[1], "r");
  if (fp == 0) {
    printf("Ooops... il file non c'è\n");
    exit(1);
  }
  for(int i =0;i<nbins; i++){
    for(int k =0;k<6; k++){
      fscanf(fp, "%lf",&dato);
      a[i][k]=dato;
	  //cout<<a[i][k]<<" "<<endl;
    }
 
  
    flux[i]=a[i][2];
    //cout<<flux[i]<<endl;
    fluxerr[i]=a[i][3];
    //    flag[i]=a[i][5];
    tstart[i]=a[i][0];
    tstop[i]=a[i][1];
    sp[i]=2.21;
    //    sperr[i]=a[i][8];
    //  cout<<sp[i]<<" "<<ra[i]<<endl;
 
 aaa=0;
 //b=1.e-18;
 ee=1000;
 freq=2.42e22*ee/100;
 emax=100000;
 emin=100;

 ums[i] = 1.-sp[i];
 conv[i]=ums[i]*pow(ee,(-sp[i]))/(pow(emax,ums[i])-pow(emin,ums[i]))*6.62e-2*(ee/100);
 // conv is in Jy
 // now convert in nufnu erg/cm2/s
 convjy[i]=conv[i]*freq*1.e-23;
 nufnu[i]=flux[i]*convjy[i];
 b[i]=fluxerr[i]*convjy[i];
 mjdsta[i]=((tstart[i]- 220838400 )/86400 )+54466;
 mjdsto[i]=((tstop[i]- 220838400 )/86400 )+54466;
 // err_log=log10((nufnu[i]+b[i])/nufnu[i]);
 

   cout<<freq<<" | "<<aaa<<" | "<<nufnu[i]<<" | "<<b[i]<<" | "<<mjdsta[i]<<" | "<<mjdsto[i]<<" |  |" <<endl;
 
 // cout<<log10(freq)<<" "<<log10(a)<<" "<<log10(nufnu)<<" "<<err_log<<endl;
  }
}
