#include "beta.h"

using namespace std;

int main(int argc, char *argv[]){

  std::string fileName;

  int tt,index=0,i1=0,ipr=0,imr=0,ipz=0,imz=0;

  int save;
 
  int pS = 5;

  static float data[5000000000],num;

  static float near[50000000];

  static float nearXi[10000000];

  static float nearT[10000000];

  float x,y,z,xi,t,r;

  float meanXi=0,squareXi=0,meanT=0,squareT=0,meanXipr=0,squareXipr=0,meanXimr=0,squareXimr=0,meanXipz=0,squareXipz=0,meanXimz=0,squareXimz=0;

  float rt = 0.9, zt = 1.5, dr = 0.00625, dz = 0.00625;

  tt = atoi(argv[1]);

  save = atoi(argv[2]);

  fileName = "output_" + std::to_string(tt) + ".dat";

  freopen(fileName.c_str(),"r", stdin);

  while(scanf("%a",&num) == 1) {
    data[index]=num;
    index++;
  }

  for ( int ii=0; ii<index/pS; ii++){
    x = data[0+pS*ii];
    y = data[1+pS*ii];
    z = data[2+pS*ii];
    xi = data[3+pS*ii];
    t = data[4+pS*ii];
    r = sqrt(x*x + y*y);

    if ( (abs(r-rt) < dr) & (abs(z-zt) < dz) ){
      near[0+pS*i1]=x;
      near[1+pS*i1]=y;
      near[2+pS*i1]=z;
      near[3+pS*i1]=xi;
      near[4+pS*i1]=t;
      i1++;
    }

  }

  if ( save == 1){
    tt++;
    fileName = "output_" + std::to_string(tt) + ".dat";
    std::ofstream myfile;
    myfile.open (fileName);
    
    for ( int ii=0;ii<i1;ii++ ){
      myfile << std::to_string( near[0 + pS*ii] ) << " " <<
                std::to_string( near[1 + pS*ii] ) << " " <<
                std::to_string( near[2 + pS*ii] ) << " " <<
                std::to_string( near[3 + pS*ii] ) << " " <<
                std::to_string( near[4 + pS*ii] ) << "\n";
    }
    myfile.close();

  }
  else {
    int ipr=0,imr=0,ipz=0,imz=0;
    for ( int ii=0;ii<i1;ii++){
      xi = near[3+pS*ii];
      nearXi[ii] = xi;
      meanXi += xi;
      squareXi += xi*xi;
      t = near[4+pS*ii];
      nearT[ii] = t;
      meanT += t;
      squareT += t*t;
      r = sqrt(near[0+pS*ii]*near[0+pS*ii] + near[1+pS*ii]*near[1+pS*ii]);

      if ( r > rt ){
        meanXipr += xi;
        squareXipr += xi*xi;
        ipr++;}
      else {
        meanXimr += xi;
        squareXimr += xi*xi;
        imr++;}
             
      if ( near[2+pS*ii] > zt ){
        meanXipz += xi;
        squareXipz += xi*xi;
        ipz++;}
      else {
        meanXimz += xi;
        squareXimz += xi*xi;
        imz++;}    
      }

      meanXi = meanXi/float(i1);
      squareXi = squareXi/float(i1);
      meanT = meanT/float(i1);
      squareT = squareT/float(i1);
      meanXipr = meanXipr/float(ipr);
      squareXipr = squareXipr/float(ipr);
      meanXimr = meanXimr/float(imr);
      squareXimr = squareXimr/float(imr);
      meanXipz = meanXipz/float(ipz);
      squareXipz = squareXipz/float(ipz);
      meanXimz = meanXimz/float(imz);
      squareXimz = squareXimz/float(imz);

      std::cout << i1 << std::endl;

      std::vector<float> xiVector(nearXi, nearXi+i1-1);

      std::sort(xiVector.begin(),xiVector.end());

      std::vector<float> tVector(nearT, nearT+i1-1);

      std::sort(tVector.begin(),tVector.end());

      fileName = "results.dat";
      std::ofstream myfile;
      myfile.open (fileName);

      myfile << std::to_string( meanXi ) << "\n";
      myfile << std::to_string( squareXi ) << "\n";
      myfile << std::to_string( meanXipr ) << "\n";
      myfile << std::to_string( squareXipr ) << "\n";
      myfile << std::to_string( meanXimr ) << "\n";
      myfile << std::to_string( squareXimr ) << "\n";
      myfile << std::to_string( meanXipz ) << "\n";
      myfile << std::to_string( squareXipz ) << "\n";
      myfile << std::to_string( meanXimz ) << "\n";
      myfile << std::to_string( squareXimz ) << "\n";
      myfile << std::to_string( rt ) << "\n";
      myfile << std::to_string( dr ) << "\n";
      myfile << std::to_string( zt ) << "\n";
      myfile << std::to_string( dz ) << "\n";

      for (std::vector<float>::iterator it=xiVector.begin(); it!=xiVector.end(); ++it){
          myfile << std::to_string( *it ) << "\n";
        }

      for (std::vector<float>::iterator it=tVector.begin(); it!=tVector.end(); ++it){
          myfile << std::to_string( *it ) << "\n";
        }

      int nx=200;
      int nt = 200;
      double Q[nx];

      //PSD_1D(nx,meanXi,squareXi,Q);
      
      PSD_Gamma(nx,nt,meanXi,squareXi,meanT,squareT,Q);

      for ( int ii=0; ii<nx; ii++ ){ myfile << std::to_string( Q[ii] ) << "\n"; }

      beta(nx,meanXi,squareXi,Q);

      for ( int ii=0; ii<nx; ii++ ){ myfile << std::to_string( Q[ii] ) << "\n"; }

      myfile.close();
  }


}


















