/****************************************************************/
/*** Main program for experiments with planar graphs on [0,1]² square     									      ***/
/*** Results (histogram in .dat files) in folder OutputPlanar_Normed	  										  ***/
/*** For evaluation and fitting with gnuplot use: load "Plotscripts/gnuplot_ausgabe_fit_powerLaw_planar.plt"      ***/
/*** For evaluation of mean shortest path use: load "Plotscripts/gnuplot_ausgabe_fit_powerLaw_planarMean.plt"     ***/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "graphs_lists.h"

int main(int argc, char **argv)
{
    int num_nodes;                            
    gs_graph_t *g;
    int num_real, i;                               
    int argz = 1;                   

	//read command line arguments
    if(argc != 2)
    {
        printf("USAGE %s <N> <m> \n", argv[0]);
        printf("N = Anzahl der Knoten im Graph\n");
        exit(1);
    }
    num_nodes = atoi(argv[argz++]);

	//run experiments in order to get distributions of shortest path in histograms for defined graph size
    runExperimentsPlanar(100,num_nodes);
	//run experiments in order to get mean of shortest path for different graph sizes
	runMeanExperimentsPlanar(10,1);
    return(0);

}
