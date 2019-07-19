#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef DEBUG
#define DEBUG_PRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( false )
#else
#define DEBUG_PRINT(...) do{ } while ( false )
#endif

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

template<class T>
void upmin(T &x, const T &y){if(x > y) x = y;}

void fprintarr(FILE * stream, int **arr, int n,int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++)
            fprintf(stream, "%10d ", arr[i][j]);
        fprintf(stream, "\n");
    }
}

int main(int argc, char **argv){
    int n_rows, n_cols;
    int **M, **E, **F;
    char *infilename = (char*)"result/input.txt";
    char *stdfilename = (char*)"result/std.txt";
    FILE *pInfile, *pStdfile;

    assert(argc == 3);
    n_rows = atoi(argv[1]);
    n_cols = atoi(argv[2]);
    M = alloc2d<int>(n_rows, n_cols);
    E = alloc2d<int>(n_rows, n_cols);
    F = alloc2d<int>(n_rows, n_cols);

    pInfile = fopen(infilename, "w");
    pStdfile = fopen(stdfilename, "w");
    
    for(int i = 0; i < n_rows; i++)
        for(int j = 0; j < n_cols; j++)
            M[i][j] = rand()*10 + rand();
 
    fprintarr(pInfile, M, n_rows, n_cols);

    for(int i = 0; i < n_rows; i++)
        for(int  j = 0; j < n_cols; j++){
            int tmp = 0;
            if (j + 1 < n_cols)
                tmp += abs(M[i][j] - M[i][j+1]);
            if (i + 1 < n_rows)
                tmp += abs(M[i][j]- M[i+1][j]);
            if (j + 1 < n_cols && i + 1 < n_rows)
                tmp += abs(M[i][j] - M[i+1][j+1]);
            E[i][j] = tmp;
        }
    
#ifdef DEBUG
    fprintarr(stderr, E, n_rows, n_cols);
#endif // DEBUG
    

    for(int j = 0; j < n_cols; j++)
        F[n_rows - 1][j] = E[n_rows - 1][j]; 
    for(int i = n_rows - 2; i >= 0; i--)
        for(int j = 0; j < n_cols; j++){
            F[i][j] = F[i+1][j];
            if(j -1 >= 0) upmin(F[i][j], F[i+1][j-1]);
            if(j+1 < n_cols) upmin(F[i][j], F[i+1][j+1]);
        }
    
    fprintarr(pStdfile, F, n_rows, n_cols);

    free2d(M);
    free2d(E);
    free2d(F);

    return 0;
}
