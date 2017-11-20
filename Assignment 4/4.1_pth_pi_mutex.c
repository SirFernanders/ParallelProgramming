#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "timer.h"

long thread_count;
long long n;
double sum;
pthread_mutex_t mutex;

void* Thread_sum(void* rank);

/* Only executed by main thread */
double Serial_pi(long long n);

int main(int argc, char* argv[]) {
    long       thread;  /* Use long in case of a 64-bit system */
    pthread_t* thread_handles;
    double start, finish;

    /* Get number of threads from command line */
    thread_count = strtol(argv[1], NULL, 10);
    n = strtoll(argv[2], NULL, 10);

    thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t));
    pthread_mutex_init(&mutex, NULL);
    sum = 0.0;


    GET_TIME(start);//Get Start Time

    for (thread = 0; thread < thread_count; thread++) {
        pthread_create(&thread_handles[thread], NULL, Thread_sum, (void *) thread);
    }

    for (thread = 0; thread < thread_count; thread++) {
        pthread_join(thread_handles[thread], NULL);
    }

    GET_TIME(finish);//Get End Time



    sum = 4.0*sum;

    printf("\nWith n = %lld terms and %ld threads,\n\n", n, thread);
    printf("The of estimate of is pi = %.15f\n", sum);
    printf("The elapsed time is %e seconds\n", finish - start);
    GET_TIME(start);
    sum = Serial_pi(n);
    GET_TIME(finish);
    printf("\nFor Single Thread\nThe estimate is pi  = %.15f\n", sum);
    printf("The elapsed time is %e seconds\n", finish - start);

    pthread_mutex_destroy(&mutex);
    free(thread_handles);

    return 0;
}  /* main */

/*------------------------------------------------------------------*/
void* Thread_sum(void* rank) {
    long my_rank = (long) rank;
    double factor;
    long long i;
    long long my_n = n/thread_count;
    long long my_first_i = my_n*my_rank;
    long long my_last_i = my_first_i + my_n;
    double my_sum = 0.0;

    if (my_first_i % 2 == 0)
        factor = 1.0;
    else
        factor = -1.0;

    for (i = my_first_i; i < my_last_i; i++, factor = -factor) {
        my_sum += factor/(2*i+1);
    }
    pthread_mutex_lock(&mutex);
    sum += my_sum;
    pthread_mutex_unlock(&mutex);

    return NULL;
}  /* Thread_sum */

double Serial_pi(long long n) {
    double sum = 0.0;
    long long i;
    double factor = 1.0;

    for (i = 0; i < n; i++, factor = -factor) {
        sum += factor/(2*i+1);
    }
    return 4.0*sum;

}  /* Serial_pi */



