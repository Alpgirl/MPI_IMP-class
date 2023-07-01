
#ifndef _MPI_H
#define _MPI_H

#include<vector>
#include <math.h>
#include "_vector4.h"


// Class for mpi implementation

class MPI_IMP {

public:
	/*MPI_IMP(int _argc, char** _argv) : argc(_argc), argv(_argv) {}
	int argc;
	char** argv;
	int size, rank;
	void init(int _argc, char** _argv);*/
	MPI_IMP(int _rank, int _size) : rank(_rank), size(_size) {}
	int rank, size;
	int N_x, N_y;
	int N_start, N_end;
	
	void print();
	void getNy(int, int);
	void messagePass(double*, int, int);
	double reduceTimeStep(double, int);
	void barrier();
	void gather(double*, double*, int, int);
	double timeNow();
};


#endif
