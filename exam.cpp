#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

char **alloc2d(const int n, const int m);
void free2d(char **p);
void printarr(const char **const arr, const int n, const int m, const char *pref);
void eachprintarr(const char **const arr, const int n, const int m, const int myrank, const int size);

const int datasize = 2;
const int halosize = 1;

void distributeBySend(const char **const global, const int globalrows, const int globalcols,
                      const int localrows, const int localcols,
                      const int rank, const int size,
                      MPI_Comm cartcomm, const int dims[2], const int coords[2]) {

    MPI_Request reqs[dims[0]*dims[1]];
    const int tag = 1;

    if (rank == 0) {
        MPI_Datatype block;
        int starts[2] = {0,0};
        int subsizes[2] = {localrows, localcols};
        int sizes[2] = {globalrows, globalcols};
        MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &block);
        MPI_Type_commit(&block);

        int reqno=0;
        for (int i=0; i<dims[0]; i++) {
            int startrow = i*datasize;
            int destcoords[2];
            destcoords[0] = i;
            for (int j=0; j<dims[1]; j++) {
                int startcol = j*datasize;
                destcoords[1] = j;

                int dest;
                MPI_Cart_rank(cartcomm, destcoords, &dest);
                MPI_Isend(&(global[startrow][startcol]), 1, block, dest, tag, cartcomm, &reqs[reqno++]);
            }
        }
    }

    char **local = alloc2d(localrows, localcols);
    MPI_Recv(&(local[0][0]), localrows*localcols, MPI_CHAR, 0, tag, cartcomm, MPI_STATUS_IGNORE);

    if (rank == 0)
        MPI_Waitall(dims[0]*dims[1], reqs, MPI_STATUS_IGNORE);

    eachprintarr((const char **const)local, localrows, localcols, rank, size);
}

void scatterAndExchange(const char **const global, const int globalrows, const int globalcols,
                      const int localrows, const int localcols,
                      const int rank, const int size,
                      MPI_Comm cartcomm, const int dims[2], const int coords[2]) {

    const int lefttag=1, righttag=2, uptag=3, downtag=4;

    char **local = alloc2d(localrows, localcols);
    for (int i=0; i<localrows; i++)
        for (int j=0; j<localcols; j++)
            local[i][j] = '.';

    MPI_Datatype tmp, globalblock;
    MPI_Datatype localblock;

    /* send just the interior data to the processors */
    int starts[2] = {0,0};
    int subsizes[2] = {datasize, datasize};
    int sizes[2] = {globalrows, globalcols};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &tmp);
    MPI_Type_create_resized(tmp, 0, sizeof(char), &globalblock);
    MPI_Type_commit(&globalblock);

    starts[0] = halosize; starts[1] = halosize;
    sizes[0] = localrows; sizes[1] = localcols;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &localblock);
    MPI_Type_commit(&localblock);

    int displs[size];
    int counts[size];
    for (int dest=0; dest<size; dest++) {
        int destcoords[2];
        MPI_Cart_coords(cartcomm, dest, 2, destcoords);
        int row = halosize + destcoords[0]*datasize;
        int col = halosize + destcoords[1]*datasize;

        counts[dest] = 1;
        displs[dest] = col + row*globalcols;
    }

    const char *ptr = (rank == 0 ? &(global[0][0]) : NULL ) ;

    MPI_Scatterv( ptr, counts, displs, globalblock,
                 &(local[0][0]), 1, localblock,
                 0, cartcomm);

    /* now send overlap data to neighbors above and below */
    int left, right, up, down;
    MPI_Cart_shift(cartcomm, 1, 1, &left, &right);
    MPI_Cart_shift(cartcomm, 0, 1, &down, &up);

    MPI_Sendrecv(&(local[localrows-2*halosize][0]), halosize*localcols, MPI_CHAR, up, uptag,
                 &(local[0][0]),                    halosize*localcols, MPI_CHAR, down, uptag,
                 cartcomm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&(local[halosize][0]),           halosize*localcols, MPI_CHAR, down, downtag,
                 &(local[localrows-halosize][0]), halosize*localcols, MPI_CHAR, up, downtag,
                 cartcomm, MPI_STATUS_IGNORE);

    /* now send overlap data to neighbors left and right */
    MPI_Datatype column;
    sizes[0] = localrows; sizes[1] = localcols;
    subsizes[0] = localrows; subsizes[1] = halosize;
    starts[0] = 0; starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &column);
    MPI_Type_commit(&column);

    MPI_Sendrecv(&(local[0][localcols-2*halosize]), 1, column, right, righttag,
                 &(local[0][0]),                    1, column, left,  righttag,
                 cartcomm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&(local[0][halosize]),             1, column, left,  lefttag,
                 &(local[0][localcols-halosize]),   1, column, right, lefttag,
                 cartcomm, MPI_STATUS_IGNORE);

    eachprintarr((const char **const)local, localrows, localcols, rank, size);
}

int main(int argc, char **argv) {
    int rank, size;
    int dims[2] = {0,0};
    int coords[2];
    int periods[2] = {0, 0};
    const int reorder = 1;
    MPI_Comm cartcomm;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
    MPI_Comm_rank(cartcomm, &rank);

    MPI_Cart_coords(cartcomm, rank, 2, coords);

    int globalcols = datasize*dims[0]+2*halosize;
    int globalrows = datasize*dims[0]+2*halosize;
    int localcols = datasize+2*halosize;
    int localrows = datasize+2*halosize;

    char **global = NULL;
    if (rank == 0) {
        global = alloc2d(globalrows, globalcols);

        for (int i=0; i<globalrows; i++)
            for (int j=0; j<globalcols; j++)
                global[i][j] = '.';

        char val = 'a';
        for (int i=halosize; i<globalrows-halosize; i++)
            for (int j=halosize; j<globalcols-halosize; j++) {
                global[i][j] = val;
                val++;
                if (val > 'z') val = 'a';
            }

        printf("Global array: ---\n");
        printarr((const char ** const)global, globalrows, globalcols, "");
    }

    if (argv[1] && !strcmp(argv[1],"sendrecv")) {
        if (rank == 0)
            printf("---\nDistributing with Send/Recv:---\n");

        distributeBySend((const char **const) global, globalrows, globalcols,
                          localrows, localcols,
                          rank, size,
                          cartcomm, dims, coords);
    } else {
        if (rank == 0)
            printf("---\nDistributing with Scatter/exchange:---\n");

        scatterAndExchange((const char **const)global, globalrows, globalcols,
                          localrows, localcols,
                          rank, size,
                          cartcomm, dims, coords);
    }

    MPI_Finalize();
    return 0;
}

char **alloc2d(const int n, const int m) {
    char *data = malloc( n*m * sizeof(int) );
    char **ptrs = malloc( n*sizeof(int *) );
    for (int i=0; i<n; i++)
        ptrs[i] = &(data[i*m]);

    return ptrs;
}

void free2d(char **p) {
    free(p[0]);
    free(p);
}

void printarr(const char **const arr, const int n, const int m, const char *pref) {
    for (int i=0; i<n; i++) {
        printf("%s", pref);
        for (int j=0; j<m; j++) printf("%c", arr[i][j]);
        printf("\n");
    }
}

void eachprintarr(const char **const arr, const int n, const int m, const int myrank, const int size) {
    char line[12];
    sprintf(line, "Rank %3d: ", myrank);
    for (int rank=0; rank<size; rank++) {
        if (rank == myrank) {
            printf("---\n");
            printarr(arr, n, m, line);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}