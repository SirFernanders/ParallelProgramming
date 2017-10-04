//
// Created by Fernando Montes on 10/3/17.
//

/* File:      histogram.c
 * Purpose:   Build a histogram from some random data
 *
 * Compile:   gcc -g -Wall -o histogram histogram.c
 * Run:       ./histogram <bin_count> <min_meas> <max_meas> <data_count>
 *
 * Input:     None
 * Output:    A histogram with X's showing the number of measurements
 *            in each bin
 *
 * Notes:
 * 1.  Actual measurements y are in the range min_meas <= y < max_meas
 * 2.  bin_counts[i] stores the number of measurements x in the range
 * 3.  bin_maxes[i-1] <= x < bin_maxes[i] (bin_maxes[-1] = min_meas)
 * 4.  DEBUG compile flag gives verbose output
 * 5.  The program will terminate if either the number of command line
 *     arguments is incorrect or if the search for a bin for a
 *     measurement fails.
 *
 * IPP:  Section 2.7.1 (pp. 66 and ff.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Usage(char prog_name[]);

void Get_args(
        char*    argv[]        /* in  */,
        int*     bin_count_p   /* out */,
        float*   min_meas_p    /* out */,
        float*   max_meas_p    /* out */,
        int*     data_count_p  /* out */,
        int*     local_data_count_p,
        int      my_rank,
        int      comm_sz,
        MPI_Comm comm
);

void Gen_data(
        float   min_meas    /* in  */,
        float   max_meas    /* in  */,
        float   data[]      /* out */,
        int     data_count  /* in  */,
        float   local_data[],
        int     local_data_count,
        int my_rank,
        MPI_Comm comm
);


void Gen_bins(
        float min_meas      /* in  */,
        float max_meas      /* in  */,
        float bin_maxes[]   /* out */,
        int   local_bin_count[]  /* out */,
        int   bin_count     /* in  */);

int Which_bin(
        float    data         /* in */,
        float    bin_maxes[]  /* in */,
        int      bin_count    /* in */,
        float    min_meas     /* in */);

void Print_histo(
        float    bin_maxes[]   /* in */,
        int      bin_counts[]  /* in */,
        int      bin_count     /* in */,
        float    min_meas      /* in */);

int main(int argc, char* argv[]) {
    int bin_count, i, bin;
    float min_meas, max_meas;
    float* bin_maxes;
    int* bin_counts;
    int* local_bin_count;
    int data_count;
    int local_d_count;
    float* data;
    float* local_d;
    int my_rank;
    int comm_sz;
    MPI_Comm comm;

    MPI_Init(NULL,NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    /* Check and get command line args */
    if (argc != 5) Usage(argv[0]);
    Get_args(argv, &bin_count, &min_meas, &max_meas, &data_count, &local_d_count, my_rank, comm_sz, comm);

    /* Allocate arrays needed */
    bin_maxes = malloc(bin_count*sizeof(float));
    bin_counts = malloc(bin_count*sizeof(int));
    local_bin_count = malloc(bin_count*sizeof(int));
    local_d = malloc(bin_count*sizeof(float));
    data = malloc(data_count*sizeof(float));

    /* Generate the data */
    Gen_data(min_meas, max_meas, data, data_count,local_d, local_d_count, my_rank, comm);

    /* Create bins for storing counts */
    Gen_bins(min_meas, max_meas, bin_maxes, local_bin_count, bin_count);

    /* Count number of values in each bin */
    for (i = 0; i < local_d_count; i++) {
        bin = Which_bin(local_d[i], bin_maxes, bin_count, min_meas);
        local_bin_count[bin]++;
    }

#  ifdef DEBUG
    printf("bin_counts = ");
   for (i = 0; i < bin_count; i++)
      printf("%d ", bin_counts[i]);
   printf("\n");
#  endif

    MPI_Reduce(local_bin_count,bin_counts,bin_count,MPI_INT,MPI_SUM,0,comm);

    /* Print the histogram */
    if(my_rank==0)
        Print_histo(bin_maxes, bin_counts, bin_count, min_meas);

    free(local_bin_count);
    free(local_d);
    free(data);
    free(bin_maxes);
    free(bin_counts);
    MPI_Finalize();
    return 0;

}  /* main */


/*---------------------------------------------------------------------
 * Function:  Usage
 * Purpose:   Print a message showing how to run program and quit
 * In arg:    prog_name:  the name of the program from the command line
 */
void Usage(char prog_name[] /* in */) {
    fprintf(stderr, "usage: %s ", prog_name);
    fprintf(stderr, "<bin_count> <min_meas> <max_meas> <data_count>\n");
    exit(0);
}  /* Usage */


/*---------------------------------------------------------------------
 * Function:  Get_args
 * Purpose:   Get the command line arguments
 * In arg:    argv:  strings from command line
 * Out args:  bin_count_p:   number of bins
 *            min_meas_p:    minimum measurement
 *            max_meas_p:    maximum measurement
 *            data_count_p:  number of measurements
 */
void Get_args(
        char*    argv[]        /* in  */,
        int*     bin_count_p   /* out */,
        float*   min_meas_p    /* out */,
        float*   max_meas_p    /* out */,
        int*     data_count_p  /* out */,
        int*     local_data_count_p,
        int      my_rank,
        int      comm_sz,
        MPI_Comm comm
) {

    if (my_rank == 0) {
        *bin_count_p = strtol(argv[1], NULL, 10);
        *min_meas_p = strtof(argv[2], NULL);
        *max_meas_p = strtof(argv[3], NULL);
        *data_count_p = strtol(argv[4], NULL, 10);


#  ifdef DEBUG
        printf("bin_count = %d\n", *bin_count_p);
               printf("min_meas = %f, max_meas = %f\n", *min_meas_p, *max_meas_p);
               printf("data_count = %d\n", *data_count_p);
#  endif
    }
    *local_data_count_p = *data_count_p / comm_sz;
    *data_count_p = *local_data_count_p * comm_sz;

    MPI_Bcast(bin_count_p,1,MPI_INT,0,comm);
    MPI_Bcast(min_meas_p,1,MPI_FLOAT,0,comm);
    MPI_Bcast(max_meas_p,1,MPI_FLOAT,0,comm);
    MPI_Bcast(data_count_p,1,MPI_INT,0,comm);
    MPI_Bcast(local_data_count_p,1,MPI_INT,0,comm);
}  /* Get_args */


/*---------------------------------------------------------------------
 * Function:  Gen_data
 * Purpose:   Generate random floats in the range min_meas <= x < max_meas
 * In args:   min_meas:    the minimum possible value for the data
 *            max_meas:    the maximum possible value for the data
 *            data_count:  the number of measurements
 * Out arg:   data:        the actual measurements
 */
void Gen_data(
        float   min_meas    /* in  */,
        float   max_meas    /* in  */,
        float   data[]      /* out */,
        int     data_count  /* in  */,
        float local_d[], int local_d_count,
        int my_rank,
        MPI_Comm comm) {
    int i;

    if (my_rank==0) {
        srandom(0);
        for (i = 0; i < data_count; i++)
            data[i] = min_meas + (max_meas - min_meas) * random() / ((double) RAND_MAX);

#  ifdef DEBUG
        printf("data = ");
               for (i = 0; i < data_count; i++)
                  printf("%4.3f ", data[i]);
               printf("\n");
#  endif
    }

    MPI_Scatter(data, local_d_count,MPI_FLOAT,local_d, local_d_count,MPI_FLOAT,0,comm);
    if(my_rank==0)
        free(data);
}  /* Gen_data */


/*---------------------------------------------------------------------
 * Function:  Gen_bins
 * Purpose:   Compute max value for each bin, and store 0 as the
 *            number of values in each bin
 * In args:   min_meas:   the minimum possible measurement
 *            max_meas:   the maximum possible measurement
 *            bin_count:  the number of bins
 * Out args:  bin_maxes:  the maximum possible value for each bin
 *            bin_counts: the number of data values in each bin
 */
void Gen_bins(
        float min_meas      /* in  */,
        float max_meas      /* in  */,
        float bin_maxes[]   /* out */,
        int   bin_counts[]  /* out */,
        int   bin_count     /* in  */) {
    float bin_width;
    int   i;

    bin_width = (max_meas - min_meas)/bin_count;

    for (i = 0; i < bin_count; i++) {
        bin_maxes[i] = min_meas + (i+1)*bin_width;
        bin_counts[i] = 0;
    }

#  ifdef DEBUG
    printf("bin_maxes = ");
   for (i = 0; i < bin_count; i++)
      printf("%4.3f ", bin_maxes[i]);
   printf("\n");
#  endif
}  /* Gen_bins */


/*---------------------------------------------------------------------
 * Function:  Which_bin
 * Purpose:   Use binary search to determine which bin a measurement
 *            belongs to
 * In args:   data:       the current measurement
 *            bin_maxes:  list of max bin values
 *            bin_count:  number of bins
 *            min_meas:   the minimum possible measurement
 * Return:    the number of the bin to which data belongs
 * Notes:
 * 1.  The bin to which data belongs satisfies
 *
 *            bin_maxes[i-1] <= data < bin_maxes[i]
 *
 *     where, bin_maxes[-1] = min_meas
 * 2.  If the search fails, the function prints a message and exits
 */
int Which_bin(
        float   data          /* in */,
        float   bin_maxes[]   /* in */,
        int     bin_count     /* in */,
        float   min_meas      /* in */) {
    int bottom = 0, top =  bin_count-1;
    int mid;
    float bin_max, bin_min;

    while (bottom <= top) {
        mid = (bottom + top)/2;
        bin_max = bin_maxes[mid];
        bin_min = (mid == 0) ? min_meas: bin_maxes[mid-1];
        if (data >= bin_max)
            bottom = mid+1;
        else if (data < bin_min)
            top = mid-1;
        else
            return mid;
    }

    /* Whoops! */
    fprintf(stderr, "Data = %f doesn't belong to a bin!\n", data);
    fprintf(stderr, "Quitting\n");
    exit(-1);
}  /* Which_bin */


/*---------------------------------------------------------------------
 * Function:  Print_histo
 * Purpose:   Print a histogram.  The number of elements in each
 *            bin is shown by an array of X's.
 * In args:   bin_maxes:   the max value for each bin
 *            bin_counts:  the number of elements in each bin
 *            bin_count:   the number of bins
 *            min_meas:    the minimum possible measurment
 */
void Print_histo(
        float  bin_maxes[]   /* in */,
        int    bin_counts[]  /* in */,
        int    bin_count     /* in */,
        float  min_meas      /* in */) {
    int i, j;
    float bin_max, bin_min;

    for (i = 0; i < bin_count; i++) {
        bin_max = bin_maxes[i];
        bin_min = (i == 0) ? min_meas: bin_maxes[i-1];
        printf("%.3f-%.3f:\t", bin_min, bin_max);
        for (j = 0; j < bin_counts[i]; j++)
            printf("X");
        printf("\n");
    }
}  /* Print_histo */

//
// Created by Fernando Montes on 10/2/17.
//

#include <stdio.h>

/* We'll be using MPI routines, definitions, etc. */
#include <mpi.h>

/* Calculate local integral  */
double Trap(double left_endpt, double right_endpt, int trap_count,
            double base_len);

/* Function we're integrating */
double f(double x);

void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p,
               int* n_p);

int main(void) {
    int my_rank, comm_sz, n, local_n;
    double a, b, h, local_a, local_b;
    double local_int, total_int;


    /* Let the system do what it needs to start up MPI */
    MPI_Init(NULL, NULL);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    Get_input(my_rank, comm_sz, &a, &b, &n);


    h = (b-a)/n;          /* h is the same for all processes */
    local_n = n/comm_sz;  /* So is the number of trapezoids  */

    /* Length of each process' interval of
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    local_int = Trap(local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f = %.15e\n",
               a, b, total_int);
    }

    /* Shut down MPI */
    MPI_Finalize();

    return 0;
} /*  main  */


/*------------------------------------------------------------------
 * Function:     Trap
 * Purpose:      Serial function for estimating a definite integral
 *               using the trapezoidal rule
 * Input args:   left_endpt
 *               right_endpt
 *               trap_count
 *               base_len
 * Return val:   Trapezoidal rule estimate of integral from
 *               left_endpt to right_endpt using trap_count
 *               trapezoids
 */
double Trap(
        double left_endpt  /* in */,
        double right_endpt /* in */,
        int    trap_count  /* in */,
        double base_len    /* in */) {
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt))/2.0;
    for (i = 1; i <= trap_count-1; i++) {
        x = left_endpt + i*base_len;
        estimate += f(x);
    }
    estimate = estimate*base_len;

    return estimate;
} /*  Trap  */


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
    return x*x;
} /* f */


void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p,
               int* n_p) {
    int dest;

    if (my_rank == 0) {
        printf("Enter a, b, and n \nExample 0 3 4000\n");
        scanf("%lf %lf %d", a_p, b_p, n_p);
        for (dest = 1; dest < comm_sz; dest++) {
            MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        }
    } else { /* my_rank != 0 */
        MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
}
