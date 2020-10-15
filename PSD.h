#include <iostream>
#include <cstdlib>
#include <random>
#include <stdio.h>
#include <math.h>
//#include <mpi.h>
#include <algorithm>
#include <vector>
#include <bits/stdc++.h>

void PSD_1D( int nx, double mean, double square, double * Q ){

  double squareQ,dQ[nx],fac=0.25;

  for ( int ii=0; ii<nx; ii++ ){

    if ( double(ii+1) < mean*double(nx) ){
      Q[ii] = 1.0;}
    else if ( double(ii) < mean*double(nx) ){
      Q[ii] = mean*double(nx) - floor(mean*double(nx));}
    else {Q[ii]=0.0;}
  }  // initialization

  squareQ = 0.0;
  for ( int ii=0; ii<nx; ii++){squareQ+=Q[ii]*Q[ii];}
  squareQ = squareQ/double(nx);

  while (squareQ > square){

    std::cout << mean << " " << squareQ << " " << square << std::endl;
    for ( int ii=0; ii<nx; ii++ ){dQ[ii]=0.0;}
    for ( int ii=1; ii<nx; ii++ ){dQ[ii]+=fac*(Q[ii-1]-Q[ii]);}
    for ( int ii=0; ii<nx-1; ii++ ){dQ[ii]+=fac*(Q[ii+1]-Q[ii]);}
    for ( int ii=0; ii<nx; ii++ ){Q[ii]+=dQ[ii];}
  
    squareQ = 0.0;
    for ( int ii=0; ii<nx; ii++){squareQ+=Q[ii]*Q[ii];}
    squareQ = squareQ/double(nx);
  }
}

bool cmp(const std::vector<double> &a,const std::vector<double> &b) 
{ 
	return a[0]<b[0]; 
}


void PSD_Gamma( int nx, int nt, double mean, double square, double meanT, double squareT, double * Q ){

  double squareQ,dQ[nx*nt],Qt[nx*nt],fac=0.25;

  double facG[nt],facT[nt],sumG,t,k,theta,mean1,var1;

  for ( int jj=0; jj<nt; jj++ ){

  for ( int ii=0; ii<nx; ii++ ){

    if ( double(ii+1) < mean*double(nx) ){
      Qt[ii+jj*nx] = 1.0;}
    else if ( double(ii) < mean*double(nx) ){
      Qt[ii+jj*nx] = mean*double(nx) - floor(mean*double(nx));}
    else {Qt[ii]=0.0;}
  }  // initialization
  }

  var1 = (squareT/(meanT*meanT) - 1.0);
  mean1 = 1.0;

  theta = var1/mean1;
  k = mean1/theta;

  std::cout << "parameters: "<< var1 << " " << mean1 << " " << k << " " << theta << std::endl;

  for ( int jj = 0; jj<nt; jj++){

    t = 0.25 + 3.0*double(jj+1)/double(nt);
  
    facT[jj] = pow(t/3.0,1.0);
    facG[jj] = pow(t,k-1.0)*exp(-t/theta);

  }

  squareQ = 0.0; sumG = 0.0;
  for ( int jj=0; jj<nt; jj++){
    for ( int ii=0; ii<nx; ii++){squareQ+=facG[jj]*Qt[ii+jj*nx]*Qt[ii+jj*nx]; sumG+=facG[jj];}  
  }
  squareQ = squareQ/sumG;

  while (squareQ > square){

    for ( int jj=0; jj<nt; jj++){
      for ( int ii=0; ii<nx; ii++ ){dQ[ii+jj*nx]=0.0;}
      for ( int ii=1; ii<nx; ii++ ){dQ[ii+jj*nx]+=facT[jj]*fac*(Qt[ii-1+jj*nx]-Qt[ii+jj*nx]);}
      //dQ[0+jj*nx]+=facT[jj]*fac*(Q[0]-Q[1]);
      for ( int ii=0; ii<nx-1; ii++ ){dQ[ii+jj*nx]+=facT[jj]*fac*(Qt[ii+1+jj*nx]-Qt[ii+jj*nx]);}
      //dQ[nx-1+jj*nx]+=facT[jj]*fac*(Q[nx-1]-Q[nx-2]);
      for ( int ii=0; ii<nx; ii++ ){Qt[ii+jj*nx]+=dQ[ii+jj*nx];}
    }
  
    squareQ = 0.0;
    for ( int jj=0; jj<nt; jj++){
      for ( int ii=0; ii<nx; ii++){squareQ+=facG[jj]*Qt[ii+jj*nx]*Qt[ii+jj*nx];}
    }
    squareQ = squareQ/sumG;
  
  std::cout << squareQ << " " << square << std::endl;}

  std::vector<std::vector<double>> QQ;

  for ( int jj=0; jj<nt; jj++){
    for ( int ii=0; ii<nx; ii++){
      QQ.push_back({Qt[ii+jj*nx],facG[jj]});
  }}

  std::sort(QQ.begin(),QQ.end(),cmp);

  double runSum=0.0;

  for ( int ii=0; ii<nx*nt; ii++ ){
    int jj = int(floor((runSum/sumG)*double(nx) - 0.5));
    runSum+=QQ[ii][1];
    int kk = int(floor((runSum/sumG)*double(nx) - 0.5));

    if (kk!=jj){
    //if (ii%nt==0){
      Q[kk] = QQ[ii][0];
      //Q[ii/nt] = QQ[ii][0];
      }    

  }
  

}












