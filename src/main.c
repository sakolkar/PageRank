#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include "timer.h"
#include "mpi.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85
#define THRESHOLD 0.0001

int main (int argc, char* argv[]){
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre, *r_local;
    int i, j;
    double damp_const;
    int iterationcount = 0;
    double start, end;

    int comm_sz, my_rank;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    int lower_bound = nodecount*my_rank / comm_sz;
    int upper_bound = nodecount*(my_rank+1) / comm_sz;
    int nodecount_local = nodecount / comm_sz;

    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    r_local = malloc((nodecount_local) * sizeof(double));

    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;


    GET_TIME(start);
    // CORE CALCULATION
    do{ 

       ++iterationcount;
        vec_cp(r, r_pre, nodecount);

        for ( i = lower_bound; i < upper_bound; ++i){
            r_local[i-lower_bound] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r_local[i-lower_bound] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            r_local[i-lower_bound] *= DAMPING_FACTOR;
            r_local[i-lower_bound] += damp_const;
        }

        MPI_Allgather(r_local, nodecount_local, MPI_DOUBLE, r, nodecount_local, MPI_DOUBLE, comm);
        
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);
    GET_TIME(end);

    MPI_Finalize();
    // post processing
    if(my_rank == 0) {
        Lab4_saveoutput(r, nodecount, end-start);
    }

    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);

    return 0;
}
