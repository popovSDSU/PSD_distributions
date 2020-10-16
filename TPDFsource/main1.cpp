// particle Lagrangian code

#include "Routines.H"

int main(int argc, char *argv[]){

std::string fileName;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
MPI_Comm_rank(MPI_COMM_WORLD, &myProc);

if ( numProcs != npx*npy*npz ){
  std::cout << "Wrong number of processors for domain decomposition: " << numProcs << "!=" << npx << "*" << npy << "*" << npz << std::endl;
  exit(1);
  }

mypx = myProc%npx; mypy = (myProc/npx)%npy; mypz = (myProc/(npx*npy))%npz;

initializePtcls();

for ( int tt=0; tt < 3001; tt++){

  evolveSDE();

  passParticles();

  particleUpkeep();

  applyMixing();
  
  if (myProc==0){std::cout << tt << " " << lastUsed << " " << firstUnused << std::endl;}

  if ( tt%100 == 0 ){
    fileName = "output_" + std::to_string(tt) + ".dat";
    outputStats(fileName);}

  }

  MPI_Finalize();

}
