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

int main(int argc, char **argv)
{
    int mode;
    int num_nodes,max_nodes,step_size;
    double k_0;
    gs_graph_t *g;
    int num_real, m, i;
    int argz = 1;

    //read command-line arguments
    if(argc <  2)
    {
        printf("Please specify experiment mode for %s:\n" , argv[0]);
        printf("mode: 1-Scale Free 2-Scale Free Constant Pref 3-Planar 4-Planar Mean\n");
        exit(1);
    }
    mode = atoi(argv[argz++]);

    if (mode == 1)
    {
        //run experiments for fixed graph size and defined parameter with normal pref attachment (path distribution)
        if (argc >= 4)
        {
            num_nodes = atoi(argv[argz++]);
            m = atoi(argv[argz++]);
            runExperiments(1000,num_nodes,m);
        }
        else
        {
            printf("USAGE %s <mode> <N> <m> \n", argv[0]);
            printf("mode = 1-Scale Free\n");
            printf("N = Anzahl der Knoten im Graph\n");
            printf("m = Anzahl der Kanten die hinzugefügt werden\n");
        }
    }
    else if (mode == 2)
    {
        //run experiments for fixed graph size and defined parameter with constant k_0 in pref attachment (path distribution)
        if (argc >= 5)
        {
            num_nodes = atoi(argv[argz++]);
            m = atoi(argv[argz++]);
            k_0 = atof(argv[argz++]);
            runExperimentsConstant(1000,num_nodes,m,k_0);
        }
        else
        {
            printf("USAGE %s <mode> <N> <m> <k0>\n", argv[0]);
            printf("mode = 2-Scale Free Constant Pref\n");
            printf("N = Anzahl der Knoten im Graph\n");
            printf("m = Anzahl der Kanten die hinzugefügt werden\n");
            printf("k_0 = Konstante für preferential attachment\n");
        }
    }
    else if (mode == 3)
    {
        //run experiments on planar graph for fixed graph size (path distribution)
        if (argc >=3)
        {
            num_nodes = atoi(argv[argz++]);
            runExperimentsPlanar(1000,num_nodes);
        }
        else
        {
            printf("USAGE %s <mode> <N>\n", argv[0]);
            printf("mode = 3-Planar Graph \n");
            printf("N = Anzahl der Knoten im Graph\n");
        }
    }
    else if (mode == 4)
    {
        //run experiments on planar graph for different graph sizes (mean shortest path length)
        if (argc >=4)
        {
            max_nodes = atoi(argv[argz++]);
            step_size = atoi(argv[argz++]);
            runMeanExperimentsPlanar(max_nodes,step_size);
        }
        else
        {
            printf("USAGE %s <mode> <max> <stepSize>\n", argv[0]);
            printf("mode = 4-Planar Graph Mean\n");
            printf("max = maximale Graphgroesse\n");
            printf("stepSize = Schrittweite für Graphgroesse\n");
        }
    }
    return(0);
}
