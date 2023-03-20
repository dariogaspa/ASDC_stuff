#include<iostream>
#include<fstream>
#include<string>
using namespace std;

int main(int argc, char *argv[]) {
  
  FILE *fp;
  FILE *fp1;
  double dato;
  double dato1;
  double lr_vero[20000];
  double lr_fake[20000];
  double lr_val_vero[20000];
  double lr_val_fake[20000];
  double rel[20000];
  double a[20000][2];
  double a1[20000][2];
  int nbins;
  double cfact[2000];
  char nomefile[100];
  ofstream out("reliability.qdp");
  fp = fopen(argv[1], "r");
  if (fp == 0) {
    printf("Ooops... il file non c'è");
    exit(1);
  }
  
fp1 = fopen(argv[2], "r");
  if (fp1 == 0) {
    printf("Ooops... 2 il file non c'è");
    exit(1);
  }

 nbins=atoi(argv[3]);

for(int i =0;i<nbins; i++){
    for(int k =0;k<2; k++){
      fscanf(fp, "%lf",&dato);
      a[i][k]=dato;
    }

    lr_vero[i]=a[i][0];
    lr_val_vero[i]=a[i][1];
    
 }

for(int u =0;u<nbins; u++){
    for(int kk =0;kk<2; kk++){
      fscanf(fp1, "%lf",&dato1);
      a1[u][kk]=dato1;
    }

    lr_fake[u]=a1[u][0];
    lr_val_fake[u]=a1[u][1];
    
 }

for (int z=0; z<nbins; z++){
  cout<<lr_vero[z]<<" "<<lr_val_vero[z]<<" "<<lr_val_fake[z]<<endl;
  //rel[z]=1-(lr_val_fake[z]/(lr_val_vero[z]+lr_val_fake[z]));
	rel[z]=1-(lr_val_fake[z]/lr_val_vero[z]);
    out<<lr_vero[z]<<"  "<<rel[z]<<endl;
  
 }
 fclose(fp1);
 fclose(fp);
 
}
