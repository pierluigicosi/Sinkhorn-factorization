/*   MyMPI.h
 *
 *   Header file for a library of matrix/vector
 *   input/output/redistribution functions.
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 4 September 2002
 */

/************************* MACROS **************************/

#define DATA_MSG           0
#define PROMPT_MSG         1
#define RESPONSE_MSG       2

#define OPEN_FILE_ERROR    -1
#define MALLOC_ERROR       -2
#define TYPE_ERROR         -3

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define PTR_SIZE           (sizeof(void*))
#define CEILING(i,j)       (((i)+(j)-1)/(j))

/***************** MISCELLANEOUS FUNCTIONS *****************/

void *my_malloc (int id, double bytes);
void  terminate (int, char *);

/*************** DATA DISTRIBUTION FUNCTIONS ***************/
void divide_array (int, int, int, int**, int**);
void divide_matrix (int id, int p, int n, int **count, int **disp);
