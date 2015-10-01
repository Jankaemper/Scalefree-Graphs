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

void exportGraphDot(gs_graph_t *g)
{
	int i =0;
	FILE *file;
	file = fopen("Output/graph_visu.dot", "w");
	fprintf(file, "graph scale_free {\n");

	for(i=0; i<g->num_nodes; i++)
	{
		elem_t *next_neighb;
		next_neighb = g->node[i].neighbors;
		while(next_neighb != NULL) {
			fprintf(file,"%d -- %d\n", i, next_neighb->info);
			next_neighb = next_neighb->next;
		}
	}
	fprintf(file, "}");
}

void printDistances(double **dist, int size)
{
	int i,j =0;
	FILE *file;
	file = fopen("Output/distance.dat", "w");

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			fprintf(file,"%d, ", (int)dist[i][j]);
		}
		fprintf(file, "\n");
	}
}

void computeShortestPaths(gs_graph_t *g,double **dist)
{
	int i,j,k =0;
	//set diagonal (loops have dist 0)
	for (i = 0; i<g->num_nodes;i++)
	{

		for (j = 0; j<g->num_nodes;j++)
		{

			if (i==j)
			{
				dist[i][j] = 0;
			}
			else
			{
				dist[i][j] = 10000000 ;
			}	
		}
	}
	//set path distances for edges (all 1)
	for(i=0; i<g->num_nodes; i++)
	{
		elem_t *next_neighb;
		next_neighb = g->node[i].neighbors;
		while(next_neighb != NULL) 
		{
			dist[i][next_neighb->info] = 1.0;	      
			next_neighb = next_neighb->next;
		}
	}
	//start actual floid warshall after inits
	for ( k=0; k<g->num_nodes;k++)
	{
		for ( i=0; i<g->num_nodes;i++)
		{
			for ( j=0; j<g->num_nodes;j++)
			{
				if (dist[i][j] > dist[i][k]  + dist[k][j])
				{
					dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
	}
}



int main(int argc, char **argv)
{
	int num_nodes;                             /* number of nodes in graph */
	gs_graph_t *g;
	int num_real, m, i;                                /* number of realizations */
	int argz = 1;                   /* for treating command line arguments */

	if(argc != 3)
	{
		printf("USAGE %s <N> <m> \n", argv[0]);
		exit(1);
	}

	num_nodes = atoi(argv[argz++]);
	m = atoi(argv[argz++]);               /* read number of nodes */
	//  sscanf(argv[argz++], "%lf", &c2);     /* read inter-group connectivity */

	g = gs_create_graph(num_nodes);
	gs_preferential_attachment(g, m);

	exportGraphDot(g);

	double **dist; 
	dist = (double**)malloc(sizeof(double*)*(g->num_nodes));
	for (i=0;i<g->num_nodes;i++)
	{

		dist[i] = (double*)malloc(sizeof(double)*g->num_nodes);
	}
	computeShortestPaths(g,dist);
	printDistances(dist,g->num_nodes);

	return(0);

}


