#include "Data.H"
#include <iostream>
#include <fstream>

double IC( double x, double y, double z ){

  double r = sqrt(x*x + y*y);

  double ans = 1.0;

  if ( r > 0.4 ){ans = 0.0;}

  return ans;

}

double bump(double x, double x1, double x2, double x3, double x4){

  double ans;

  if ( x < x1 ){ans=0.0;}
  else if ( x < x2 ){ans = sin(0.5*M_PI*(x-x1)/(x2-x1))*sin(0.5*M_PI*(x-x1)/(x2-x1));}
  else if ( x < x3 ){ans = 1.0;}
  else if ( x < x4 ){ans = cos(0.5*M_PI*(x-x3)/(x4-x3))*cos(0.5*M_PI*(x-x3)/(x4-x3));}
  else { ans = 0.0;}

  return ans;

}

void getU(double x, double y, double z, double *U){


  double z1 = z + 2, r = pow(x*x + y*y + z1*z1, 1.5);

  U[0] = 2*x/r; U[1] = 2*y/r; U[2] = 2*z1/r;  

  return;
  

}

void getGamma(double x, double y, double z, double *Gam){

  double r = sqrt(x*x + y*y);

  *Gam = 0.015*bump(r,-0.95,-0.85,1.05,1.15)*bump(z,0.05,0.15,2.25,2.35);

  return;

}

void Drift( double x, double y, double z, double *D ){

  double dG[3], Gam1, Gam2, dd = 0.0001;

  getU(x,y,z,D);

  getGamma(x+dd,y,z,&Gam1);getGamma(x-dd,y,z,&Gam2); D[0]+=(Gam1 - Gam2)/(2.0*dd);
  getGamma(x,y+dd,z,&Gam1);getGamma(x,y-dd,z,&Gam2); D[1]+=(Gam1 - Gam2)/(2.0*dd);
  getGamma(x,y,z+dd,&Gam1);getGamma(x,y,z-dd,&Gam2); D[2]+=(Gam1 - Gam2)/(2.0*dd);

  return;

}

void Diffusion( double x, double y, double z, double *G ){

  getGamma(x,y,z,G);

  *G = sqrt(2.0*(*G));

  return;

}

void findPreviousUsed( int *lastUsed1 ){

  while ( (used[*lastUsed1] == 0) ){

  *lastUsed1 = *lastUsed1-1;

  if ( *lastUsed1 < 0 ){

    std::cout << "Error in findPreviousUsed" << std::endl;
    exit(1);

  }}

}

void findNextUnused( int *firstUnused1 ){

  while ( (used[*firstUnused1] == 1) ){
  
  *firstUnused1 = *firstUnused1+1;

  if ( *firstUnused1 > maxPtcl ){

  std::cout << "Error in findNextUnused " <<  *firstUnused1 << " " << maxPtcl << std::endl;
  exit(1);

  }}

}

void initializePtcls(){

  double x,y,z;

  for (int kk=0;kk<nz;kk++){for (int jj=0;jj<ny;jj++){for (int ii=0;ii<nx;ii++){

  for ( int pp=0;pp<ppc;pp++ ){

    lastUsed++; firstUnused++;
    x = randU(gen)*dx + double(mypx*nx + ii - (npx*nx)/2)*dx;
    y = randU(gen)*dy + double(mypy*ny + jj - (npy*ny)/2)*dy;
    z = randU(gen)*dz + double(mypz*nz + kk)*dz;

    ptcl[0 + lastUsed*pS]=x; ptcl[1 + lastUsed*pS]=y; ptcl[2 + lastUsed*pS]=z;

    ptcl[3 + lastUsed*pS] = IC(x,y,z);

    ptcl[4 + lastUsed*pS] = 0.0;

    used[lastUsed] = 1;

  }

  }}}


}

void evolveSDE(){

  int iPtcl = 0;
  double D[3],dW[3],G,x,y,z;

  while ( iPtcl <= lastUsed ){

    if ( used[iPtcl] == 1){

      x = ptcl[0 + iPtcl*pS]; y = ptcl[1 + iPtcl*pS]; z = ptcl[2 + iPtcl*pS];

      Drift(x,y,z,D); Diffusion(x,y,z,&G);

      dW[0] = sqrtDt*randN(gen);
      dW[1] = sqrtDt*randN(gen);
      dW[2] = sqrtDt*randN(gen);

      x += D[0]*dt + G*dW[0];
      y += D[1]*dt + G*dW[1];
      z += D[2]*dt + G*dW[2];

      ptcl[0 + iPtcl*pS] = x; ptcl[1 + iPtcl*pS] = y; ptcl[2 + iPtcl*pS] = z;

      ptcl[4 + iPtcl*pS] +=dt;

    }

    iPtcl++;

  }

}


void addParticlesFromBuf( int howMany ){

  for ( int ii=0; ii<howMany; ii++){
    ptcl[0 + firstUnused*pS] = recvBuffer[0+ii*pS];
    ptcl[1 + firstUnused*pS] = recvBuffer[1+ii*pS];
    ptcl[2 + firstUnused*pS] = recvBuffer[2+ii*pS];
    ptcl[3 + firstUnused*pS] = recvBuffer[3+ii*pS];
    ptcl[4 + firstUnused*pS] = recvBuffer[4+ii*pS];

    used[firstUnused] = 1;

    if ( firstUnused > lastUsed ){lastUsed = firstUnused;}


    findNextUnused(&firstUnused);

  }

}


void passParticles(){

  int howMany[1];
  int iPtcl;
  double x,y,z;
  MPI_Status Stat;

// pass in the positive x-direction
  if (mypx > 0){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc-1,0,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc-1,1,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypx < npx-1){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      x = ptcl[0 + iPtcl*pS];
      if ( (x >= double(mypx*nx+nx-(nx*npx)/2)*dx) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc+1,0,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc+1,1,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the positive x-direction

// pass in the negative x-direction
  if (mypx < npx-1){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc+1,2,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc+1,3,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypx > 0){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      x = ptcl[0 + iPtcl*pS];
      if ( (x < double(mypx*nx-(nx*npx)/2)*dx) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc-1,2,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc-1,3,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the negative x-direction

// pass in the positive y-direction
  if (mypy > 0){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc-npx,4,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc-npx,5,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypy < npy-1){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      y = ptcl[1 + iPtcl*pS];
      if ( (y >= double(mypy*ny+ny-(ny*npy)/2)*dy) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc+npx,4,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc+npx,5,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the positive y-direction

// pass in the negative y-direction
  if (mypy < npy-1){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc+npx,6,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc+npx,7,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypy > 0){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      y = ptcl[1 + iPtcl*pS];
      if ( (y < double(mypy*ny-(ny*npy)/2)*dy) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc-npx,6,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc-npx,7,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the negative y-direction

// pass in the positive z-direction
  if (mypz > 0){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc-npx*npy,8,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc-npx*npy,9,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypz < npz-1){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      z = ptcl[2 + iPtcl*pS];
      if ( (z >= double(mypz*nz+nz)*dz) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc+npx*npy,8,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc+npx*npy,9,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the positive z-direction

// pass in the negative z-direction
  if (mypz < npz-1){
    MPI_Recv(howMany,1,MPI_INTEGER,myProc+npx*npy,10,MPI_COMM_WORLD,&Stat);
    if ( *howMany > 0 ){
      MPI_Recv(recvBuffer,pS*(*howMany),MPI_DOUBLE,myProc+npx*npy,11,MPI_COMM_WORLD,&Stat);
      addParticlesFromBuf(*howMany);
    }
  } // recv particles first

  if (mypz > 0){
    *howMany = 0; iPtcl=0;
    while (iPtcl <= lastUsed){if (used[iPtcl] == 1){
      z = ptcl[2 + iPtcl*pS];
      if ( (z < double(mypz*nz)*dz) ){
        used[iPtcl]=0;
        sendBuffer[pS*(*howMany)+0] = ptcl[0 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+1] = ptcl[1 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+2] = ptcl[2 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+3] = ptcl[3 + iPtcl*pS];
        sendBuffer[pS*(*howMany)+4] = ptcl[4 + iPtcl*pS];
        *howMany = *howMany + 1;
      }      
    } iPtcl++; }
    MPI_Send(howMany,1,MPI_INTEGER,myProc-npx*npy,10,MPI_COMM_WORLD);

    if ( *howMany > 0 ){    
      MPI_Send(sendBuffer,pS*(*howMany),MPI_DOUBLE,myProc-npx*npy,11,MPI_COMM_WORLD);
      findPreviousUsed(&lastUsed);
      firstUnused=0;
      findNextUnused(&firstUnused);
      }
  } // send particles next
// end pass in the negative z-direction

}


void particleUpkeep(){

  double x, y, z;
  int iPtcl = 0;

  // disable out of bounds ptcls
  
  while (iPtcl <= lastUsed){

    if (used[iPtcl] == 1){
      x = ptcl[0 + iPtcl*pS]; y = ptcl[1 + iPtcl*pS]; z = ptcl[2 + iPtcl*pS];

      if ( (z < double(mypz*nz)*dz)|(z >= double(mypz*nz + nz)*dz)|(x < double(mypx*nx-(nx*npx)/2)*dx)|(x >= double(mypx*nx+nx-(nx*npx)/2)*dx)|
           (y < double(mypy*ny-(ny*npy)/2)*dy)|(y >= double(mypy*ny + ny-(ny*npy)/2)*dy) ){ 
            used[iPtcl] = 0; 
            //std::cout << "particle disabled" << std::endl;
            }
      if ((mypz==0)&(z<dz)){used[iPtcl] = 0;} // reseed particles near z=0
    }

    iPtcl++; 

    }

  findPreviousUsed(&lastUsed);
  firstUnused=0;
  findNextUnused(&firstUnused);  

  // end disable out of bounds ptcls
  

  // introduce new ptcls
  
  if (mypz==0){
  for (int jj=0;jj<ny;jj++){for (int ii=0;ii<nx;ii++){

  for ( int pp=0;pp<ppc;pp++ ){

    x = randU(gen)*dx + double(mypx*nx + ii - (nx*npx)/2)*dx;
    y = randU(gen)*dy + double(mypy*ny + jj - (ny*npy)/2)*dy;
    z = randU(gen)*dz;

    ptcl[0 + firstUnused*pS]=x; ptcl[1 + firstUnused*pS]=y; ptcl[2 + firstUnused*pS]=z;

    ptcl[3 + firstUnused*pS] = IC(x,y,z);
    
    ptcl[4 + firstUnused*pS] = 0.0;

    used[firstUnused] = 1;

    if ( firstUnused > lastUsed ){lastUsed = firstUnused;}

    findNextUnused(&firstUnused);

  }

  }}}


  // end introduce new ptcls

  
  // sort the ptcls
  while ( lastUsed > firstUnused ){

    used[lastUsed] = 0; used[firstUnused] = 1;
    for ( int ii=0;ii<pS;ii++ ){ ptcl[ii + pS*firstUnused] = ptcl[ii + pS*lastUsed];}

    findPreviousUsed(&lastUsed);
    findNextUnused(&firstUnused);

    //std::cout << lastUsed << " " << firstUnused << std::endl;

  }
  // end sort the ptcls

}

bool cmp(const std::vector<double> &a,const std::vector<double> &b) 
{ 
	return a[1]<b[1]; 
}

void mixBin(int count, double * bPtcl, double Omega){

  /*/ IEM
  
  double mmean=0;

    for ( int ii=0; ii<count; ii++){mmean += bPtcl[ii];}

    mmean = mmean/double(count);

    for ( int ii=0; ii<count; ii++){

      bPtcl[ii] = mmean + (bPtcl[ii]-mmean)*pow(Euler,-Omega*dt);

    }

  // end IEM */

  // EMST
  
  double dphi[count],dx,Fo;

  int idx,times;
  
  std::vector<std::vector<double>> D;
  

  for ( int ii=0; ii<count; ii++ ){

    D.push_back({double(ii) + 0.5, bPtcl[ii]});

  }

  std::sort(D.begin(),D.end(),cmp); 

  Fo = Omega*dt*double(count)*double(count);

  times=1;

  if ( Fo > 0.25 ){times=int(ceil(Fo/0.25) + 0.1);}

  for ( int kk=0; kk<times; kk++){

    for ( int ii=0;ii<count; ii++ ){
      dphi[ii]=0;

      if ( ii > 0 ){dphi[ii]+=(D[ii-1][1] - D[ii][1]);}
      if ( ii < count-1 ){dphi[ii]+=(D[ii+1][1] - D[ii][1]);}

      dphi[ii]=Omega*dt*dphi[ii]*(double(count)*double(count)/double(times));

      }

    for ( int ii=0;ii<count; ii++ ){D[ii][1]+=dphi[ii];}

  }

  for ( int ii=0; ii<count; ii++ ){
    idx = int(D[ii][0]);
    bPtcl[idx] = D[ii][1]; 
  }

  // end EMST

}

void applyMixing(){

  int i1,j1,k1,idx1,maxCount=0,minCount=2*pS*ppc;
  double Omega = 0.1;

  for (int kk=0;kk<nz;kk++){for (int jj=0;jj<ny;jj++){for (int ii=0;ii<nx;ii++){
    binCount[ii + nx*jj + nx*ny*kk] = 0;
  }}}

  std::cout << "qqq1" << std::endl;

  for (int ii=0; ii<=lastUsed; ii++){
    i1 = int(floor(ptcl[ 0 + pS*ii ]/dx)) - mypx*nx + (nx*npx)/2;
    j1 = int(floor(ptcl[ 1 + pS*ii ]/dy)) - mypy*ny + (ny*npy)/2;
    k1 = int(floor(ptcl[ 2 + pS*ii ]/dz)) - mypz*nz;

    if ( (i1 < 0)|(i1 >= nx) ){
	std::cout << "Bad i1: " << i1 << " " << nx;
        exit(1);}

    if ( (j1 < 0)|(j1 >= ny) ){
	std::cout << "Bad j1: " << j1 << " " << ny;
        exit(1);}

    if ( (k1 < 0)|(k1 >= nz) ){
	std::cout << "Bad k1: " << k1 << " " << nz;
        exit(1);}

    idx1 = binCount[i1 + nx*j1 + nx*ny*k1];
    binCount[i1 + nx*j1 + nx*ny*k1] = idx1+1;
    binIdx[ idx1 + 2*pS*ppc*(i1 + nx*j1 + nx*ny*k1) ] = ii;
    binPtcl[ idx1 + 2*pS*ppc*(i1 + nx*j1 + nx*ny*k1) ] = ptcl[ 3 + pS*ii ];
  }

  for (int kk=0;kk<nz;kk++){for (int jj=0;jj<ny;jj++){for (int ii=0;ii<nx;ii++){

    if (binCount[ii+nx*jj+nx*ny*kk] > maxCount){maxCount = binCount[ii+nx*jj+nx*ny*kk];}
    if (binCount[ii+nx*jj+nx*ny*kk] < minCount){minCount = binCount[ii+nx*jj+nx*ny*kk];}

    mixBin(binCount[ii+nx*jj+nx*ny*kk],&binPtcl[0 + 2*pS*ppc*(ii+nx*jj+nx*ny*kk)], Omega);

    for ( int pp=0; pp<binCount[ii+nx*jj+nx*ny*kk]; pp++ ){
      ptcl[3 + pS*binIdx[ pp + 2*pS*ppc*(ii+nx*jj+nx*ny*kk)]] = 
          binPtcl[ pp + 2*pS*ppc*(ii+nx*jj+nx*ny*kk)];
      }

  }}}

  //std::cout << "min,max count: " << minCount << " " << maxCount << std::endl;

}

void outputStats(const std::string& fileName){

  std::ofstream myfile;

for ( int zz=0; zz<numProcs; zz++){
if ( zz==myProc ){

  if (zz==0){
    myfile.open (fileName);}
  else{
    myfile.open (fileName,std::ios_base::app);
  }

  for (int ii=0; ii <= lastUsed; ii++){

    myfile << std::to_string( ptcl[0 + pS*ii] ) << " " <<
              std::to_string( ptcl[1 + pS*ii] ) << " " <<
              std::to_string( ptcl[2 + pS*ii] ) << " " <<
              std::to_string( ptcl[3 + pS*ii] ) << " " <<
              std::to_string( ptcl[4 + pS*ii] ) << "\n"; 

    }


  myfile.close();

}

MPI_Barrier(MPI_COMM_WORLD);

}

}






















