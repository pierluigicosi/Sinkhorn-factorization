#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "MyMPI.h"


/***************** MISCELLANEOUS FUNCTIONS *****************/

/*
 *   Given MPI_Datatype 't', function 'get_size' returns the
 *   size of a single datum of that data type.
 */

int get_size (MPI_Datatype t) {
   if (t == MPI_BYTE) return sizeof(char);
   if (t == MPI_DOUBLE) return sizeof(double);
   if (t == MPI_FLOAT) return sizeof(float);
   if (t == MPI_INT) return sizeof(int);
   printf ("Error: Unrecognized argument to 'get_size'\n");
   fflush (stdout);
   MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}


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

void scatter_arrays (
   int id,       /* IN - Process rank */
   int p,        /* IN - Number of processes */
   int n,        /* IN - Total number of elements */
   int **count,  /* OUT - Array of counts */
   int **disp)   /* OUT - Array of displacements */
{

   int i;

   *count = malloc (p * sizeof(int));
   *disp = malloc (p * sizeof(int));
   (*count)[0] = BLOCK_SIZE(0,p,n)*n;
   (*disp)[0] = 0;
   for (i = 1; i < p; i++) {
      (*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      (*count)[i] = BLOCK_SIZE(i,p,n)*n;
   }
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

void create_mixed_xfer_arrays (
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

void create_uniform_xfer_arrays (
   int id,        /* IN - Process rank */
   int p,         /* IN - Number of processes */
   int n,         /* IN - Number of elements */
   int **count,   /* OUT - Array of counts */
   int **disp)    /* OUT - Array of displacements */
{

   int i;

   *count = my_malloc (id, p * sizeof(int));
   *disp = my_malloc (id, p * sizeof(int));
   (*count)[0] = BLOCK_SIZE(id,p,n);
   (*disp)[0] = 0;
   for (i = 1; i < p; i++) {
      (*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      (*count)[i] = BLOCK_SIZE(id,p,n);
   }
}

/*
 *   This function is used to transform a vector from a
 *   block distribution to a replicated distribution within a
 *   communicator.
 */

void replicate_block_vector (
   void        *ablock,  /* IN - Block-distributed vector */
   int          n,       /* IN - Elements in vector */
   void        *arep,    /* OUT - Replicated vector */
   MPI_Datatype dtype,   /* IN - Element type */
   MPI_Comm     comm)    /* IN - Communicator */
{
   int *cnt;  /* Elements contributed by each process */
   int *disp; /* Displacement in concatenated array */
   int id;    /* Process id */
   int p;     /* Processes in communicator */

   MPI_Comm_size (comm, &p);
   MPI_Comm_rank (comm, &id);
   create_mixed_xfer_arrays (id, p, n, &cnt, &disp);
   MPI_Allgatherv (ablock, cnt[id], dtype, arep, cnt,
                   disp, dtype, comm);
   free (cnt);
   free (disp);
}

/******************** OUTPUT FUNCTIONS ********************/

/*
 *   Print elements of a doubly-subscripted array.
 */

void print_submatrix (
   void       **a,       /* OUT - Doubly-subscripted array */
   MPI_Datatype dtype,   /* OUT - Type of array elements */
   int          rows,    /* OUT - Matrix rows */
   int          cols)    /* OUT - Matrix cols */
{
   int i, j;

   for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
         if (dtype == MPI_DOUBLE)
            printf ("%6.3f ", ((double **)a)[i][j]);
         else {
            if (dtype == MPI_FLOAT)
               printf ("%6.3f ", ((float **)a)[i][j]);
            else if (dtype == MPI_INT)
               printf ("%6d ", ((int **)a)[i][j]);
         }
      }
      putchar ('\n');
   }
}


/*
 *   Print elements of a singly-subscripted array.
 */

void print_subvector (
   void        *a,       /* IN - Array pointer */
   MPI_Datatype dtype,   /* IN - Array type */
   int          n)       /* IN - Array size */
{
   int i;

   for (i = 0; i < n; i++) {
      if (dtype == MPI_DOUBLE)
         printf ("%6.3f ", ((double *)a)[i]);
      else {
         if (dtype == MPI_FLOAT)
            printf ("%6.3f ", ((float *)a)[i]);
         else if (dtype == MPI_INT)
            printf ("%6d ", ((int *)a)[i]);
      }
   }
}


/*
 *   Print a matrix distributed checkerboard fashion among
 *   the processes in a communicator.
 */

void print_checkerboard_matrix (
   void       **a,            /* IN -2D matrix */
   MPI_Datatype dtype,        /* IN -Matrix element type */
   int          m,            /* IN -Matrix rows */
   int          n,            /* IN -Matrix columns */
   MPI_Comm     grid_comm)    /* IN - Communicator */
{
   void      *buffer;         /* Room to hold 1 matrix row */
   int        coords[2];      /* Grid coords of process
                                 sending elements */
   int        datum_size;     /* Bytes per matrix element */
   int        els;            /* Elements received */
   int        grid_coords[2]; /* Coords of this process */
   int        grid_id;        /* Process rank in grid */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dims of process grid */
   int        i, j, k;
   void      *laddr;          /* Where to put subrow */
   int        local_cols;     /* Matrix cols on this proc */
   int        p;              /* Number of processes */
   int        src;            /* ID of proc with subrow */
   MPI_Status status;         /* Result of receive */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coords);
   local_cols = BLOCK_SIZE(grid_coords[1],grid_size[1],n);

   if (!grid_id)
      buffer = my_malloc (grid_id, n*datum_size);

   /* For each row of the process grid */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by the process row */
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],m); j++) {

         /* Collect the matrix row on grid process 0 and
            print it */
         if (!grid_id) {
            for (k = 0; k < grid_size[1]; k++) {
               coords[1] = k;
               MPI_Cart_rank (grid_comm, coords, &src);
               els = BLOCK_SIZE(k,grid_size[1],n);
               laddr = buffer +
                  BLOCK_LOW(k,grid_size[1],n) * datum_size;
               if (src == 0) {
                  memcpy (laddr, a[j], els * datum_size);
               } else {
                  MPI_Recv(laddr, els, dtype, src, 0,
                     grid_comm, &status);
               }
            }
            print_subvector (buffer, dtype, n);
            putchar ('\n');
         } else if (grid_coords[0] == i) {
            MPI_Send (a[j], local_cols, dtype, 0, 0,
               grid_comm);
         }
      }
   }
   if (!grid_id) {
      free (buffer);
      putchar ('\n');
   }
}


/*
 *   Print a matrix that has a columnwise-block-striped data
 *   decomposition among the elements of a communicator.
 */

void print_col_striped_matrix (
   void       **a,       /* IN - 2D array */
   MPI_Datatype dtype,   /* IN - Type of matrix elements */
   int          m,       /* IN - Matrix rows */
   int          n,       /* IN - Matrix cols */
   MPI_Comm     comm)    /* IN - Communicator */
{
   MPI_Status status;     /* Result of receive */
   int        datum_size; /* Bytes per matrix element */
   void      *buffer;     /* Enough room to hold 1 row */
   int        i, j;
   int        id;         /* Process rank */
   int        p;          /* Number of processes */
   int*       rec_count;  /* Elements received per proc */
   int*       rec_disp;   /* Offset of each proc's block */

   MPI_Comm_rank (comm, &id);
   MPI_Comm_size (comm, &p);
   datum_size = get_size (dtype);
   create_mixed_xfer_arrays (id, p, n, &rec_count,&rec_disp);

   if (!id)
      buffer = my_malloc (id, n*datum_size);

   for (i = 0; i < m; i++) {
      MPI_Gatherv (a[i], BLOCK_SIZE(id,p,n), dtype, buffer,
         rec_count, rec_disp, dtype, 0, MPI_COMM_WORLD);
      if (!id) {
         print_subvector (buffer, dtype, n);
         putchar ('\n');
      }
   }
   free (rec_count);
   free (rec_disp);
   if (!id) {
      free (buffer);
      putchar ('\n');
   }
}


/*
 *   Print a matrix that is distributed in row-striped
 *   fashion among the processes in a communicator.
 */

void print_row_striped_matrix (
   void **a,            /* IN - 2D array */
   MPI_Datatype dtype,  /* IN - Matrix element type */
   int m,               /* IN - Matrix rows */
   int n,               /* IN - Matrix cols */
   MPI_Comm comm)       /* IN - Communicator */
{
   MPI_Status  status;          /* Result of receive */
   void       *bstorage;        /* Elements received from
                                   another process */
   void      **b;               /* 2D array indexing into
                                   'bstorage' */
   int         datum_size;      /* Bytes per element */
   int         i;
   int         id;              /* Process rank */
   int         local_rows;      /* This proc's rows */
   int         max_block_size;  /* Most matrix rows held by
                                   any process */
   int         prompt;          /* Dummy variable */
   int         p;               /* Number of processes */

   MPI_Comm_rank (comm, &id);
   MPI_Comm_size (comm, &p);
   local_rows = BLOCK_SIZE(id,p,m);
   if (!id) {
      print_submatrix (a, dtype, local_rows, n);
      if (p > 1) {
         datum_size = get_size (dtype);
         max_block_size = BLOCK_SIZE(p-1,p,m);
         bstorage = my_malloc (id,
            max_block_size * n * datum_size);
         b = (void **) my_malloc (id,
            max_block_size * datum_size);
         b[0] = bstorage;
         for (i = 1; i < max_block_size; i++) {
            b[i] = b[i-1] + n * datum_size;
         }
         for (i = 1; i < p; i++) {
            MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
               MPI_COMM_WORLD);
            MPI_Recv (bstorage, BLOCK_SIZE(i,p,m)*n, dtype,
               i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
            print_submatrix (b, dtype, BLOCK_SIZE(i,p,m), n);
         }
         free (b);
         free (bstorage);
      }
      putchar ('\n');
   } else {
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG,
         MPI_COMM_WORLD, &status);
      MPI_Send (*a, local_rows * n, dtype, 0, RESPONSE_MSG,
         MPI_COMM_WORLD);
   }
}


/*
 *   Print a vector that is block distributed among the
 *   processes in a communicator.
 */

void print_block_vector (
   void        *v,       /* IN - Address of vector */
   MPI_Datatype dtype,   /* IN - Vector element type */
   int          n,       /* IN - Elements in vector */
   MPI_Comm     comm)    /* IN - Communicator */
{
   int        datum_size; /* Bytes per vector element */
   int        i;
   int        prompt;     /* Dummy variable */
   MPI_Status status;     /* Result of receive */
   void       *tmp;       /* Other process's subvector */
   int        id;         /* Process rank */
   int        p;          /* Number of processes */

   MPI_Comm_size (comm, &p);
   MPI_Comm_rank (comm, &id);
   datum_size = get_size (dtype);

   if (!id) {
      print_subvector (v, dtype, BLOCK_SIZE(id,p,n));
      if (p > 1) {
         tmp = my_malloc (id,BLOCK_SIZE(p-1,p,n)*datum_size);
         for (i = 1; i < p; i++) {
            MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
               comm);
            MPI_Recv (tmp, BLOCK_SIZE(i,p,n), dtype, i,
               RESPONSE_MSG, comm, &status);
            print_subvector (tmp, dtype, BLOCK_SIZE(i,p,n));
         }
         free (tmp);
      }
      printf ("\n\n");
   } else {
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm,
         &status);
      MPI_Send (v, BLOCK_SIZE(id,p,n), dtype, 0,
         RESPONSE_MSG, comm);
   }
}


/*
 *   Print a vector that is replicated among the processes
 *   in a communicator.
 */

void print_replicated_vector (
   void        *v,      /* IN - Address of vector */
   MPI_Datatype dtype,  /* IN - Vector element type */
   int          n,      /* IN - Elements in vector */
   MPI_Comm     comm)   /* IN - Communicator */
{
   int id;              /* Process rank */

   MPI_Comm_rank (comm, &id);

   if (!id) {
      print_subvector (v, dtype, n);
      printf ("\n\n");
   }
}
