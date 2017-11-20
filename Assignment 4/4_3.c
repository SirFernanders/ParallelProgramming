#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


long thread_count;

int flag;
message_available = 0;
consumer =1;
pthread_mutex_t mutex;


void* Thread_sum(void* rank);


int main(int argc, char* argv[]) {
    long thread;
    pthread_t* thread_handles;
    thread_count = strtol(argv[1], NULL, 10);

    thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t));
    pthread_mutex_init(&mutex,NULL);


    flag = 0;
    for (thread = 0; thread < thread_count; thread++) {
        printf("%ld \n", thread);
        pthread_create(&thread_handles[thread], NULL, Thread_sum, (void *) thread);
    }

    for (thread = 0; thread < thread_count; thread++)
        pthread_join(thread_handles[thread], NULL);


    free(thread_handles);
    return 0;
}

void* Thread_sum(void* rank) {
    
    int i;


    i=0;
    while(1){
        pthread_mutex_lock(&mutex);
        if (message_available){
            printf("thread %li : message consumed\n",rank);
            message_available = 0;
            pthread_mutex_unlock(&mutex);
        }else{
            printf("thread %li : message produced\n",rank);
            message_available =1;
            pthread_mutex_unlock(&mutex);
            break;
        }
        i++;
        pthread_mutex_unlock(&mutex);
    }



    return NULL;
}


