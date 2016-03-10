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
    double *r, *r_pre;
    int i, j;
    double damp_const;
    int iterationcount = 0;
    double start, end;

    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;

    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    GET_TIME(start);
    // CORE CALCULATION
    do{
        ++iterationcount;
        vec_cp(r, r_pre, nodecount);
        for ( i = 0; i < nodecount; ++i){
            r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            r[i] *= DAMPING_FACTOR;
            r[i] += damp_const;
        }
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);
    GET_TIME(end);

    // post processing
    Lab4_saveoutput(r, nodecount, end-start);

    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);

    return 0;
}
