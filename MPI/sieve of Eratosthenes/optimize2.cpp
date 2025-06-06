// decrease the cost of communication
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define LL long long

//优化：将marked0的大小改为sqrt(n)的大小，同时所有进程均复制一份，减少了while中的通信开销
int main(int argc, char *argv[])
{
    LL count;            /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    LL first;            /* Index of first multiple */
    LL first0;
    LL global_count = 0; /* Global prime count */
    LL high_value;       /* Highest value on this proc */
    LL i;
    int id;        /* Process ID number */
    LL index;      /* Index of current prime */
    LL low_value;  /* Lowest value on this proc */
    char *marked;  /* Portion of 2,...,'n' */
    char *marked0; /* Portion of 2,...,'sqrt(n)' */
    LL n;          /* Sieving from 2, ..., 'n' */
    int p;         /* Number of processes */
    LL proc0_size; /* Size of proc 0's subarray */
    LL prime;      /* Current prime */
    LL size;       /* Elements in 'marked' */
    LL size0;      /* Elements in 'marked0'*/
    LL low0;
    LL high0;
    LL start_n = (int)sqrt((double)atoll(argv[1])) + 1; //所有进程共享的变量为需要处理总数的开方
    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2)
    {
        if (!id)
            printf("Command line: %s <m>\n", argv[0]); 
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]); //

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */
    //最低数为除开前n开方后的偏移量，同样也是奇数索引，无偶数
    low_value = start_n + id * (n - start_n + 1) / p; // lowest value in the process
    if (low_value % 2 == 0)
    {
        low_value += 1;
    }                                                            // make sure the lowest value is odd, to calculate the right size
    high_value = start_n - 1 + (id + 1) * (n - start_n + 1) / p; // highest value in the process
    size = (high_value - low_value) / 2 + 1; // 按道理应该要high_value - low_value + 1再除2，但这里要取上限，所以-1除2再+1
    // 根号n的后一个数归在各自进程后面算，故所有共享的只需要算前根号n包括前面一个数的部分，故是-4而不是-3   
    size0 = (start_n - 4) / 2 + 1; // the range we need to find the primes used in marked size0同理

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    /* Allocate this process's share of the array. */

    marked = (char *)malloc(size);   // the array start_n~n to find the primes
    marked0 = (char *)malloc(size0); // the array 3~sqrt(n) to find var prime

    if (marked == NULL || marked0 == NULL)
    { // allocate the space for arraies
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size0; i++)
    {
        marked0[i] = 0;
        marked[i] = 0;
    } // size0 is less than size in this situation
    for (; i < size; i++)
        marked[i] = 0;

    index = 0; // initial the var
    prime = 3;
    do
    {
        // find the first multiple in marked0
        first0 = (prime * prime - 3) / 2;  //进程0中的初始位置单独找
        // find the first multiple in marked
        if (!(low_value % prime))
            first = 0;
        else
            first = prime - (low_value % prime);
        //下面if判断为了转换为奇数索引，找的都是在本进程中第一个primer的奇数倍数
        if ((low_value + first) % 2 != 0)
            first = first / 2;
        else
            first = (first + prime) / 2;
        // marked the multiples in marked0 and marked
        for (i = first0; i < size0; i += prime)
            marked0[i] = 1;
        for (i = first; i < size; i += prime)
            marked[i] = 1;
        // find the next prime
        while (marked0[++index] == 0);
        prime = index * 2 + 3;
    } while (prime * prime <= high_value);
    //前根号n的统计工作也平分给所有进程
    count = 0;
    low0 = id * size0 / p;
    high0 = (id + 1) * size0 / p - 1;
    // count the primes in marked0
    for (i = low0; i <= high0; i++)
        if (!marked0[i]){
            printf("prime %lld\n", i * 2 + 3);
            count++;
        }
            
    // count the primes in marked
    for (i = 0; i < size; i++)
        if (!marked[i]){
            printf("prime %lld\n", (i * 2 + low_value));
            count++;
        }
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
               0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id)
    {
        printf("There are %lld primes less than or equal to %lld\n",
               global_count + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}