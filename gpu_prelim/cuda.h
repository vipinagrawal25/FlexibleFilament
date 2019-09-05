#ifndef FILE_CUDA_SEEN
#define FILE_CUDA_SEEN
/*--------------------------------------------------*/
#define lmessage 2048 
struct CRASH{
  int lstop;
  char message[lmessage] ;
};
const int size_CRASH = sizeof(int) + lmessage*sizeof(char);
void  qdevice( int *count, cudaDeviceProp **prop ) ;
void qfree( cudaDeviceProp *prop );
__device__ void device_exception( struct CRASH *bug, char *mesg );
__device__ void scpy( char to[], char from[]);
void set_crash( CRASH *BUG, CRASH **dev_bug ) ;
void IStop( CRASH BUG );
__global__ void thread_maxima( double array[], double redux[]);
__global__ void thread_sum( double array[], double redux[]);
#endif /* !CUDA_SEEN */
