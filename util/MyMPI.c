#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "MyMPI.h"


/***************** MISCELLANEOUS FUNCTIONS *****************/


/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */

void *my_malloc (
   int id,     /* IN - Process rank */
   double bytes)  /* IN - Bytes to allocate */
{
   void *buffer;
    if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed for process %d\n", id);
      fflush (stdout);
      MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
   }
   return buffer;
}


/*
 *   Function 'terminate' is called when the program should
 *   not continue execution, due to an error condition that
 *   all of the processes are aware of. Process 0 prints the
 *   error message passed as an argument to the function.
 *
 *   All processes must invoke this function together!
 */

void terminate (
   int   id,            /* IN - Process rank */
   char *error_message) /* IN - Message to print */
{
   if (!id) {
      printf ("Error: %s\n", error_message);
      fflush (stdout);
   }
   MPI_Finalize();
   exit (-1);
}


/************ DATA DISTRIBUTION FUNCTIONS ******************/

/*
 *   This function creates the count and displacement arrays
 *   needed by scatter and gather functions, when the number
 *   of elements send/received to/from other processes
 *   varies.
 */

void divide_array (
   int id,       /* IN - Process rank */
   int p,        /* IN - Number of processes */
   int n,        /* IN - Total number of elements */
   int **count,  /* OUT - Array of counts */
   int **disp)   /* OUT - Array of displacements */
{

   int i;

   *count = my_malloc (id, p * sizeof(int));
   *disp = my_malloc (id, p * sizeof(int));
   (*count)[0] = BLOCK_SIZE(0,p,n);
   (*disp)[0] = 0;
   for (i = 1; i < p; i++) {
      (*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      (*count)[i] = BLOCK_SIZE(i,p,n);
   }
}


/*
 *   This function creates the count and displacement arrays
 *   needed in an all-to-all exchange, when a process gets
 *   the same number of elements from every other process.
 */

void divide_matrix (
   int id,       /* IN - Process rank */
   int p,        /* IN - Number of processes */
   int n,        /* IN - Total number of elements */
   int **count,  /* OUT - Array of counts */
   int **disp)   /* OUT - Array of displacements */
{

   int i;

   *count = my_malloc (id,p * sizeof(int));
   *disp = my_malloc (id,p * sizeof(int));
   (*count)[0] = BLOCK_SIZE(0,p,n)*n;
   (*disp)[0] = 0;
   for (i = 1; i < p; i++) {
      (*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      (*count)[i] = BLOCK_SIZE(i,p,n)*n;
   }
}

