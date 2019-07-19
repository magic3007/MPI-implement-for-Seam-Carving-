#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

using namespace std;

template<class T> 
T **alloc2d(const size_t n, const size_t m);

template<class T>
void free2d(T **p);

template<class T>
void upmin(T &x, const T &y){if(x > y) x = y;}

enum MPI_TAG_T{
    MPI_SCATTER_M,
    MPI_GATHER_F,
    MPI_COM_E
};

int main(int argc, char **argv){
    int n_rows, n_cols;
    FILE *pInfile;
    const char *filename;
    int comm_rank, comm_size;
    int **M, **E, **F;
    int *col_counts, *col_shifts;
    int **recvbuf, recvcount;
    int col_count, col_shift;
    int n_buf_cols;
    MPI_Request *reqs;
    int left_F, right_F;
    MPI_Request left_req, right_req, req;

    assert(argc==4);
    n_rows = atoi(argv[1]);
    n_cols = atoi(argv[2]);
    filename = argv[3];

    pInfile = fopen(filename, "r");

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    if(comm_rank == 0){

        // row first
        M = alloc2d<int>(n_rows, n_cols + 1);
        for(int i = 0; i < n_rows; i++)
            for(int j = 0; j < n_cols; j++)
                fscanf(pInfile, "%d", &M[i][j]);
        
        col_counts = new int[comm_size];
        col_shifts = new int[comm_size];
        for(int i = 1; i < comm_size; i++)
            col_counts[i] = n_cols / comm_size;
        col_counts[0] = n_cols - n_cols / comm_size * (comm_size - 1);
        
        col_shifts[0] = 0;
        for(int i = 1; i < comm_size; i++)
            col_shifts[i] = col_shifts[i-1] + col_counts[i-1];
    }

    MPI_Scatter(col_counts, comm_size, MPI_INT, &col_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(col_shifts, comm_size, MPI_INT, &col_shift, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    if (comm_rank == 0)
        fprintf(stderr, "n_rows=%d n_cols=%d comm_size=%d\n", n_rows, n_cols, comm_size);
    fprintf(stderr, "rank %d: col_count=%d col_shift=%d\n", comm_rank, col_count, col_shift);
#endif

    MPI_Datatype type;
    int sizes[2]    = {n_rows, n_cols};  /* size of global array */
    int subsizes[2] = {n_rows, n_cols / comm_size + 1};  /* size of sub-region */
    int starts[2]   = {0,0};  /* let's say we're looking at region "0",
                                 which begins at index [0,0] */
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
    MPI_Type_commit(&type);

    if (comm_rank == 0){
        n_buf_cols = n_cols + 1;
        recvbuf = M;
        reqs = new MPI_Request[comm_size];
        for(int i = 1; i < comm_size; i++)
            MPI_Isend(&(M[0][col_shifts[i]]), 1, type, i, MPI_SCATTER_M, MPI_COMM_WORLD, reqs + i);
        MPI_Waitall(comm_size - 1, reqs + 1, MPI_STATUS_IGNORE);
    }else{
        n_buf_cols = col_count + 1;
        recvbuf = alloc2d<int>(n_rows, n_buf_cols);
        MPI_Recv(recvbuf[0], n_rows * n_buf_cols, MPI_INT, 0, MPI_SCATTER_M, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(comm_rank == 0){
        E = alloc2d<int>(n_rows, col_count);
        F = alloc2d<int>(n_rows, n_cols);
    }else{
        E = alloc2d<int>(n_rows, col_count);
        F = alloc2d<int>(n_rows, col_count);
    }

    for(int i = 0; i < n_rows; i++)
        for(int j = col_shift; j < col_shift + col_count; j++){
            int tmp = 0;
            if ( j + 1 < n_cols)
                tmp += abs(recvbuf[i][j - col_shift] - recvbuf[i][j + 1 - col_shift]);
            if ( i + 1 < n_rows)
                tmp += abs(recvbuf[i][j - col_shift] - recvbuf[i + 1][j - col_shift]);
            if ( j + 1 < n_cols && i + 1 < n_rows)
                tmp += abs(recvbuf[i][j - col_shift] - recvbuf[i + 1][j + 1 - col_shift]);
            E[i][j - col_shift] = tmp;   
        }
        
    for(int i = n_rows - 1; i >= 0; i--){
        if( i == n_rows -1){
            for(int j = col_shift; j < col_shift + col_count; j++)
                F[i][j - col_shift] = E[i][j - col_shift];
        }else{
            if(comm_rank != 0)
                MPI_Recv(&left_F, 1, MPI_INT, comm_rank - 1, MPI_COM_E + i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(comm_rank != comm_size - 1)
                MPI_Recv(&right_F, 1, MPI_INT, comm_rank + 1, MPI_COM_E + i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j = col_shift; j < col_shift + col_count; j++){
                F[i][j - col_shift] = F[i + 1][j - col_shift];
                if(j != 0)
                    upmin(F[i][j - col_shift], j == col_shift ? left_F : F[i + 1][j - col_shift - 1]);
                if(j != n_cols - 1)
                    upmin(F[i][j - col_shift], j == col_shift + col_count - 1 ? right_F : F[i + 1][j - col_shift + 1]);
                F[i][j - col_shift] += E[i][j - col_shift];
            }
        }
        if(comm_rank != 0){
            if( i != n_rows - 1)
                MPI_Wait(&left_req, MPI_STATUS_IGNORE);
            MPI_Isend(&F[i][0], 1, MPI_INT, comm_rank - 1, MPI_COM_E + i, MPI_COMM_WORLD, &left_req);
        }
        if(comm_rank != comm_size - 1){
            if( i != n_rows - 1)
                MPI_Wait(&right_req, MPI_STATUS_IGNORE);
            MPI_Isend(&F[i][col_count - 1], 1, MPI_INT, comm_rank + 1, MPI_COM_E + i, MPI_COMM_WORLD, &right_req);
        }
    }
    
    if(comm_rank != 0){
        MPI_Send(F[0], n_rows * col_count, MPI_INT, 0, MPI_GATHER_F, MPI_COMM_WORLD);
    }else{
        for(int i = 1; i < comm_size; i++)
            MPI_Irecv(&(F[0][col_shifts[i]]), 1, type, i, MPI_GATHER_F, MPI_COMM_WORLD, reqs + i);
        MPI_Waitall(comm_size - 1, reqs + 1, MPI_STATUS_IGNORE);

        for(int i = 0; i < n_rows; i++){
            for(int j = 0; j < n_cols; j++)
                fprintf(stdout, "%d ", F[i][j]);
            printf("\n");
        }
        
    }

    if(comm_rank == 0){
        free2d(M);
        delete [] col_counts;
        delete [] col_shifts;
        delete [] reqs;
    }else{
        free2d(recvbuf);
    }
    free2d(E);
    free2d(F);
    
    MPI_Type_free(&type);
    MPI_Finalize();
    return 0;
}


template<class T> 
T **alloc2d(const size_t n, const size_t m){
    T *data = (T*)malloc(n * m * sizeof(T));
    T **ptrs = (T**)malloc(n * sizeof(T*));
    for (size_t i = 0; i < n; i++)
        ptrs[i] = &(data[i * m]);
    return ptrs;
}

template<class T>
void free2d(T **p){
    free(p[0]);
    free(p);
}
