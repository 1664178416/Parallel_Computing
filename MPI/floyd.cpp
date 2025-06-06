#include <mpi.h>
#include <stdio.h>
#include <cstdlib>

#define BLOCK_OWNER(index, p, n) ((p * (index + 1) - 1) / n)
#define BLOCK_LOW(id, p, n) (id * n / p)
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id + 1), p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id + 1), p, n) - BLOCK_LOW(id, p, n))
#define MIN(a, b) (a < b ? a : b)

void print(int *array, int n);
void printrow(int *array, int low, int high, int n);
int main(int argc, char *argv[])
{
    int n; // size of the matrix
    int i, j, k;
    int id;
    int p;
    int min, max;
    double time, max_time;
    int *a;
    int offset; // local index of broadcast row
    int root;   // process controlling row to be becast
    int *temp;  // holds the broadcast rows

    if (argc != 4)
    {
        printf("Please input two parameter :the second is the size of the matrix");
    }
    n = atoi(argv[1]);
    min = atoi(argv[2]);
    max = atoi(argv[3]); 
    a = (int *)malloc(n * n * sizeof(int));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                a[i * n + j] = 0;
            }
            else
            {
                a[i * n + j] = (rand() % (max - min + 1)) + min;
            }
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    if (!id)
    {
        printf("the origin matrix is:\n");
        print(a, n);
    }
    time = -MPI_Wtime();
    temp = (int *)malloc(n * sizeof(int));
    for (k = 0; k < n; k++)
    {
        //在本行的负责发，不在本行的负责收和计算
        root = BLOCK_OWNER(k, p, n);
        if (root == id)
        {
            offset = k - BLOCK_LOW(id, p, n);
            for (j = 0; j < n; j++)
            {
                temp[j] = a[offset * n + j];
            }
        }
        //广播和接收都是走同一条指令
        MPI_Bcast(temp, n, MPI_INT, root, MPI_COMM_WORLD);
        //n是行数，i也是行数，按行划分
        for (i = 0; i < BLOCK_SIZE(id, p, n); i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i * n + j] = MIN(a[i * n + j], a[i * n + k] + temp[j]);
            }
        }
    }
    free(temp);
    time += MPI_Wtime();
    int *recv_counts = (int *)malloc(p * sizeof(int)); // element number each process received
    int *displs = (int *)malloc(p * sizeof(int));      // displs in every process
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    for (int i = 0; i < p; i++)
    {
        recv_counts[i] = n * BLOCK_SIZE(i, p, n);
        displs[i] = n * BLOCK_LOW(i, p, n);
    }
    int *recv_buffer = NULL;
    if (id == 0)
    {
        recv_buffer = (int *)malloc(n * n * sizeof(int)); // buffer to save the receive data
    }
    MPI_Gatherv((a + n * BLOCK_LOW(id, p, n)), BLOCK_SIZE(id, p, n) * n, MPI_INT,
                recv_buffer, recv_counts, displs, MPI_INT,
                0, MPI_COMM_WORLD);
    if (!id)
    {
        printf("Floyd\nmatrix size: %d\n%dprocess:%6.2f seconds\n", n, p, max_time);
        printf("result\n");
        print(recv_buffer, n);
    }
    free(a);
    MPI_Finalize();
    return 0;
}
void print(int *array, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d,", array[i * n + j]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}
void printrow(int *array, int low, int high, int n)
{
    for (int i = low; i <= high; i++)
    {
        printf("row%d:  ", i);
        for (int j = 0; j < n; j++)
        {
            printf("%d,", array[i * n + j]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }
}
