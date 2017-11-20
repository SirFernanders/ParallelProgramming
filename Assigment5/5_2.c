//
// Created by Fernando Montes on 11/19/17.
//

#include<stdio.h>
#include<math.h>
#include<pthread.h>
#include <unistd.h>

int a[100], size, N;
int min, max;
int grno;
float grsize;
int histo[50];

int start=0;
pthread_mutex_t lock;
pthread_cond_t cond;

void * work()
{
    int i, temp;

    pthread_mutex_lock(&lock);
    if(start==0)
        pthread_cond_wait(&cond, &lock);
    pthread_mutex_unlock(&lock);

    for(;;)
    {
        pthread_mutex_lock(&lock);
        i=size-1;
        size--;
        pthread_mutex_unlock(&lock);

        if(size<0)
        {
            pthread_mutex_lock(&lock);
            pthread_cond_wait(&cond, &lock);
            pthread_mutex_unlock(&lock);
        }

        if(i<0)
            break;

        if(min>a[i])
            min=a[i];
        if(max<a[i])
            max=a[i];
    }

    for(;;)
    {
        pthread_mutex_lock(&lock);
        i=size-1;
        size--;
        pthread_mutex_unlock(&lock);

        if(i<0)
            break;

        temp=(float)((a[i]-min))/grsize;
        if(temp >(grno-1))
            temp= grno-1;
        histo[temp]++;
    }
}

int main()
{
    pthread_t id;
    void * status;
    int i, nthread;
    int temp1, temp2;

    pthread_mutex_init(&lock,0);
    pthread_cond_init(&cond,0);

    printf("Enter the size of the array :: ");
    scanf("%d",&size);
    N=size;

    printf("\nEnter the %d elements of array\n",N);
    for(i=0;i<N;i++)
        scanf("%d",&a[i]);
    min=max=a[0];

    printf("\nEnter the total number of groups required in histogram::");
    scanf("%d",&grno);

    for(i=0;i<grno;i++)
        histo[i]=0;

    printf("\nEnter the total number of threads :: ");
    scanf("%d",&nthread);
    for(i=0;i<nthread;i++)
    {
        if(0==pthread_create(&id,0,work,0))
            continue;
        else
            printf("\nError in creating threads");
    }

    pthread_mutex_lock(&lock);
    start=1;
    if(0!=pthread_cond_broadcast(&cond))
        printf("\nError in broadcasting");
    pthread_mutex_unlock(&lock);

    while(size>=0)
        sleep(1);

    grsize=(float)(max-min)/(float)grno;
    size=N;

    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&cond);
    pthread_mutex_unlock(&lock);


    while(size>=0)
        sleep(1);

    pthread_join(id,status);

    printf("\nThe histogram is as follows::\n");
    temp1=min;
    temp2=ceil(grsize);

    for(i=0;i<grno;i++)
    {
        if(i!=grno-1)
        {
            printf("For %d -> %d = %d\n",temp1, temp2,histo[i]);
            temp1=temp2+1;
            temp2+=ceil(grsize);
        }
        else
            printf("For %d and up values = %d ",temp1, histo[i]);
    }

    pthread_mutex_destroy(&lock);

    printf("\n");
    return 0;
}
