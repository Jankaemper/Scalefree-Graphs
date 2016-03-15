/****************************************************************/
/*** Main program for undirected graphs                       ***/
/*** A.K. Hartmann January 2008                               ***/
/*** Forgeschrittene Computersimulationen                     ***/
/*** University of Oldenburg, Germany 2008                    ***/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "graphs_lists.h"

void testMain(int num_nodes,int m)
{
    int i;
    gs_graph_t *g;
    g = gs_create_graph(num_nodes);
    gs_preferential_attachment(g, m);

    exportGraphDot(g,num_nodes,m);

    double **dist;
    dist = (double**)malloc(sizeof(double*)*(g->num_nodes));
    for (i=0;i<g->num_nodes;i++)
    {
        dist[i] = (double*)malloc(sizeof(double)*g->num_nodes);
    }
    computeShortestPaths(g,dist);
    printDistances(dist,num_nodes,m);
}


int main(int argc, char **argv)
{
    int num_nodes;                             /* number of nodes in graph */
    gs_graph_t *g;
    int num_real, i;                                /* number of realizations */
    int argz = 1;                   /* for treating command line arguments */

    if(argc != 2)
    {
        printf("USAGE %s <N> <m> \n", argv[0]);
        printf("N = Anzahl der Knoten im Graph\n");
        exit(1);
    }

    num_nodes = atoi(argv[argz++]);

    runExperimentsPlanar(10000,num_nodes);
    return(0);

}
