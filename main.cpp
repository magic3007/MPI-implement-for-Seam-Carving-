#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv){
    int n_rows, n_cols;
    FILE *pInfile;
    const char *filename;
    int comm_rank, comm_size;
    int *M, *E;
    int *row_counts, *row_shifts;
    int *sendcounts, *displs, *recvbuf, recvcount;
    int *sendbuf, sendcount;
    int r_count, r_shift;

    assert(argc==4);
    n_rows = atoi(argv[1]);
    n_cols = atoi(argv[2]);
    filename = argv[3];

    pInfile = fopen(filename, "r");

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    row_counts = new int[comm_size];
    row_shifts = new int[comm_size];
    sendcounts = new int[comm_size];
    displs = new int[comm_size];

    for(int i = 0 ; i < comm_size; i++)
        row_counts[i] = n_rows / comm_size;
    row_counts[0] = n_rows - (n_rows / comm_size) * comm_size;
    
    row_shifts[0] = 0;
    for(int i = 1; i < comm_size; i++)
        row_shifts[i] = row_shifts[i-1] + row_counts[i-1];
    
    for(int i = 0; i < comm_size; i++){
        sendcounts[i] = (i == comm_size - 1 ? row_counts[i] : row_counts[i] + 1) * n_cols;
        displs[i] = row_shifts[i] * n_cols;
    }

    if(comm_rank == 0){
        M = new int[n_rows * n_cols];
        for(int i = 0; i < n_rows; i++)
            for(int j = 0; j < n_cols; j++)
                fscanf(pInfile, "%d", &M[i * n_cols + j]);
    }

    recvcount = sendcounts[comm_rank];
    recvbuf = new int[recvcount];

    MPI_Scatterv(M, sendcounts, displs, MPI_INT, recvbuf, recvcount, MPI_INT, 0, MPI_COMM_WORLD);
    r_count = row_counts[comm_rank];
    r_shift = row_shifts[comm_rank];

    sendbuf = new int[r_count * n_cols];
    for(int i = r_shift; i < r_shift + r_count; i++)
        for(int j = 0; j < n_cols; j++){
            int tmp = 0;
            if(j + 1 < n_cols) 
                tmp += abs(recvbuf[(i-r_shift) * n_cols + j] - recvbuf[(i - r_shift) * n_cols + j +1 ]);
            if(i + 1 < n_rows)
                tmp += abs(recvbuf[(i - r_shift) * n_cols + j] - recvbuf[(i + 1 - r_shift) * n_cols + j]);
            if(i + 1 < n_rows && j + 1 < n_cols)
                tmp += abs(recvbuf[(i - r_shift) * n_cols + j] - recvbuf[(i + 1 - r_shift) * n_cols + j + 1]);
            sendbuf[(i - r_shift) * n_cols + j] = tmp;
        }
    
    
    MPI_Finalize();
    delete [] sendcounts;
    delete [] displs;
    return 0;
}