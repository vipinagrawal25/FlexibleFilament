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
#endif /* !CUDA_SEEN */
