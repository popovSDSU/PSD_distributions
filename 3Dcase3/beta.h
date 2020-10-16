#include "PSD.h"

void beta( int nx, double mean, double square, double * Q ){

  double sigma = square - mean*mean;
  double p = mean*(mean*(1.0-mean)/sigma -1.0);
  double q = (1.0-mean)*(mean*(1.0-mean)/sigma -1.0);
  double zp,zm,z,dQ;

  dQ = 1.0/(nx-1);

  Q[0] = 0.0;

  for ( int ii=1; ii < nx; ii++ ){
    zp = double(ii)*dQ; zm = double(ii-1)*dQ;
    z = (zp + zm)/2.0;

    if ( ii < nx/2 ){
      Q[ii] = Q[ii-1] + (1.0/p)*(pow(1-z,q-1.0))*(pow(zp,p) - pow(zm,p));
    }
    else{
      Q[ii] = Q[ii-1] - (1.0/q)*(pow(z,p-1.0))*(pow(1.0-zp,q) - pow(1.0-zm,q));
    }
  }

  for ( int ii=0; ii < nx; ii++){Q[ii]=Q[ii]/Q[nx-1];}

}
