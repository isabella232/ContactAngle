
#ifndef mingw
#include "mpi.h"
#else
typedef int MPI_Comm ;
typedef int MPI_Group;
typedef int MPI_Request;
typedef int MPI_Status;
typedef char MPI_BYTE;
typedef char MPI_Datatype;
typedef char MPI_Op;
#define  MPI_MIN  0
#define  MPI_MAX  1
#define  MPI_SUM  2
#define  MPI_LOR  3
#define  MPI_LAND 4

#define MPI_INT    0
#define MPI_LABEL  1
#define MPI_SCALAR 2
#define MPI_DOUBLE 2

#define MPI_Group_incl(...)
#define MPI_Comm_create(...)
#define MPI_Initialized(...)
#define MPI_Init(...)
#define MPI_Finalize(...)
#define MPI_Bsend(...) 0
#define MPI_Send(...) 0
#define MPI_Isend(...) 0
#define MPI_Allreduce(...)
#define MPI_Probe(...)
#define MPI_Recv(...) 0
#define MPI_Irecv(...) 0

#define MPI_COMM_WORLD 1
#define MPI_GROUP_NULL 0
#define MPI_COMM_NULL  0
#define MPI_Wait(...) 1
#define MPI_Test(...)
#define MPI_Waitall(...)  0
#define MPI_Abort(...)
#define MPI_Buffer_detach(...)
#define MPI_Comm_size(...)
#define MPI_Comm_rank(...) 0
#define MPI_Buffer_attach(...)
#define MPI_Comm_free(...)
#define MPI_Group_free(...)
#define MPI_Get_count(...)
#define MPI_Comm_group(...)


#endif



