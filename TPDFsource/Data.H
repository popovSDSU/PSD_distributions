#include <iostream>
#include <random>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <algorithm>
#include <vector>

const int nx = 20, ny = 20, nz = 30;

const double dx = 0.02, dy = 0.02, dz = 0.02;

const int maxPtcl = 12000000;

const double dt = 0.01;

double sqrtDt = sqrt(dt);

double Euler = 2.71828182845904;

const int pS = 5;

const int ppc = 60;

const int MPIbufSize = pS*ppc*nx*ny*nz/4;

static double sendBuffer[MPIbufSize];
static double recvBuffer[MPIbufSize];

int lastUsed = -1, firstUnused = 0;

double ptcl[pS*maxPtcl];
int used[maxPtcl];

int binCount[nx*ny*nz];
int binIdx[nx*ny*nz*2*pS*ppc];
double binPtcl[nx*ny*nz*2*pS*ppc];


int numProcs, npx = 6, npy = 6, npz = 4;
int myProc, mypx, mypy, mypz;

std::default_random_engine gen;
std::uniform_real_distribution<double> randU(0.0,1.0);
std::normal_distribution<double> randN(0.0,1.0);
