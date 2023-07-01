
#include "_mpi.h"
#include "mpi.h"

/*void MPI_IMP::init(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}*/
void MPI_IMP::print() {
	std::cout << "I am " << rank << " say hello!" << std::endl;
}


void MPI_IMP::getNy(int Nx, int Ny) {
    int i, cnt = 0;
	int NRemainCell;
	int* dims = new int[size];
    int* dims_shift = new int[size];
    Ny += 1;
	NRemainCell = Ny - Ny / size * size;

    for (i = 0; i < size; i++) {
        dims[i] = Ny / size + 2;
    }

    for (i = 0; i < NRemainCell; i++) {
        if (cnt == size - 1) cnt = 0;
        dims[cnt] += 1;
        cnt++;
    }
    dims[0] += 1;
    dims_shift[0] = 0;
    for (i = 1; i < size; i++) {
        dims_shift[i] = dims[i - 1] + dims_shift[i - 1];
    }
    N_x = Nx;
    N_y = dims[rank];
    N_start = dims_shift[rank];
    N_end = N_start + N_y - 1;
}


double MPI_IMP::reduceTimeStep(double var, int count) {
    double res_t;
    MPI_Reduce(&var, &res_t, count, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&res_t, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return res_t;
}


void MPI_IMP::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}



void MPI_IMP::messagePass(double* y, int end, int start) {
    MPI_Status status;
    int flag = 0;

    if (size > 1) {
        if (rank % 2 == 0) {
            if (rank < size - 1) MPI_Recv(&y[end], 1, MPI_DOUBLE, rank + 1, flag, MPI_COMM_WORLD, &status);
            if (rank < size - 1) MPI_Send(&y[end - 1], 1, MPI_DOUBLE, rank + 1, flag, MPI_COMM_WORLD);

            if (rank > 0) MPI_Send(&y[start + 1], 1, MPI_DOUBLE, rank - 1, flag, MPI_COMM_WORLD);
            if (rank > 0) MPI_Recv(&y[start], 1, MPI_DOUBLE, rank - 1, flag, MPI_COMM_WORLD, &status);
        }
        else if (rank % 2 != 0) {
            if (rank > 0) MPI_Send(&y[start + 1], 1, MPI_DOUBLE, rank - 1, flag, MPI_COMM_WORLD);
            if (rank > 0) MPI_Recv(&y[start], 1, MPI_DOUBLE, rank - 1, flag, MPI_COMM_WORLD, &status);

            if (rank < size - 1) MPI_Recv(&y[end], 1, MPI_DOUBLE, rank + 1, flag, MPI_COMM_WORLD, &status);
            if (rank < size - 1) MPI_Send(&y[end - 1], 1, MPI_DOUBLE, rank + 1, flag, MPI_COMM_WORLD);
        }
    }
    return;
}


void MPI_IMP::gather(double* y, double* U_temp, int jmin, int jmax) {
    int tmp = jmax - jmin;
    int* displs = new int[size];
    int* recvcounts = new int[size];
    displs[0] = 0;
    MPI_Gather(&tmp, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }
    }
    MPI_Gatherv(&y[0], jmax - jmin, MPI_DOUBLE, &U_temp[0], recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(displs);
    free(recvcounts);
}


double MPI_IMP::timeNow(){
    return MPI_Wtime();
}

