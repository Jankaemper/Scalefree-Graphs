/****************************************************************/
/*** Main program for experiments with scale-free graphs 			      		***/
/*** Results (histogram in .dat files) in folders Output/ and  Output_Normed/	***/
/*** For evaluation and fitting with gnuplot use: e.g. load "Plotscripts/gnuplot_ausgabe_fit_powerLaw_M2.plt"      ***/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "graphs_lists.h"

//test main used for debugging 
//prints graph into dot-readable format in order to visualize the graph layout with all nodes and edges
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
    int num_nodes;                             
    gs_graph_t *g;
    int num_real, m, i;                        
    int argz = 1;                 

	//read command-line arguments
    if(argc != 3)
    {
        printf("USAGE %s <N> <m> \n", argv[0]);
        printf("N = Anzahl der Knoten im Graph\n");
        printf("m = Anzahl der Kanten die hinzugefügt werden\n");
        exit(1);
    }
    num_nodes = atoi(argv[argz++]);
    m = atoi(argv[argz++]);               /* read number of nodes */

    //run experiments for fixed graph size and defined parameter
    runExperiments(10000,num_nodes,m);
    return(0);

}
