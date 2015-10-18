#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "graphs_lists.h"

/******************* gs_create_graph() ********************/
/** Creates a graph with a fixed number of nodes and    **/
/** no edges.                                           **/
/** PARAMETERS: (*)= return-paramter                    **/
/**      num_nodes: number of nodes                     **/
/** RETURNS:                                            **/
/**     pointer to created graph                        **/
/*********************************************************/
gs_graph_t *gs_create_graph(int num_nodes)
{
    gs_graph_t *g;
    int n;

    g = (gs_graph_t *) malloc(sizeof(gs_graph_t));  /* allocate */
    g->node = (gs_node_t *) malloc(num_nodes*sizeof(gs_node_t));

    g->num_nodes = num_nodes;                     /* initialise */
    for(n=0; n<num_nodes; n++)
	g->node[n].neighbors = NULL;

    return(g);
}



/******************* gs_insert_edge() ********************/
/** Insert undirected edge (from,to) into graph         **/
/** if the edge does not exists so far.                 **/
/** PARAMETERS: (*)= return-paramter                    **/
/**             g: graph                                **/
/**      from, to: id of nodes                          **/
/** RETURNS:                                            **/
/**     (nothing)                                       **/
/*********************************************************/
void gs_insert_edge(gs_graph_t *g, int from, int to)
{
    elem_t *elem1, *elem2;

    if(search_info(g->node[from].neighbors, to) != NULL)  /* edge exists? */
	return;

    elem1 = create_element(to);              /* create neighbor for 'from' */
    g->node[from].neighbors =
	insert_element(g->node[from].neighbors, elem1, NULL);
    elem2 = create_element(from);              /* create neighbor for 'to' */
    g->node[to].neighbors =
	insert_element(g->node[to].neighbors, elem2, NULL);
}


/******************* gs_edge_exists() ********************/
/** Tests whether undirected edge (from,to) in graph    **/
/** exists already.                                     **/
/** PARAMETERS: (*)= return-paramter                    **/
/**             g: graph                                **/
/**      from, to: id of nodes                          **/
/** RETURNS:                                            **/
/**     1 if edge exists, 0 else                        **/
/*********************************************************/
int gs_edge_exists(gs_graph_t *g, int from, int to)
{
    if(search_info(g->node[from].neighbors, to) != NULL)  /* edge exists? */
      return(1);
    else
      return(0);
}

/******************* gs_preferential_attachment() ********/
/** Generates graph bei "adding" one node after the     **/
/** other. For each "added" node, exactly m edges to    **/
/** already "old added" nodes are inserted randomly. The**/
/** probability that an "old added" nodes receives a new**/
/** edge is proportional to its current degree          **/
/** No self loops are allowed. No edge is allowed to    **/
/*+ appear twice!                                       **/
/** PARAMETERS: (*)= return-paramter                    **/
/**         (*) g: graph                                **/
/**             m: number of edges to be added          **/
/** RETURNS:                                            **/
/**     (nothing)                                       **/
/*********************************************************/
void gs_preferential_attachment(gs_graph_t *g, int m)
{
  int t;
  int n1, n2;
  int *pick;            /* array which holds for each edge {n1,n2} */
                        /* the numbers n1 and n2. Used for picking */
                       /* nodes proportional to its current degree */
  int num_pick;              /* number of entries in 'pick' so far */
  int max_pick;                       /* maximum number of entries */

  if(g->num_nodes < m+1)
  {
      printf("graph too small to have at least %d edges per node!\n", m);
      exit(1);
  }
  max_pick = 2*m*g->num_nodes- m*(m+1);
  pick = (int *) malloc(max_pick*sizeof(int));
  num_pick=0;
  for(n1=0; n1<m+1; n1++) /* start: complete subgraph of m+1 nodes */

    for(n2=n1+1; n2<m+1; n2++)
    {
      gs_insert_edge(g, n1, n2);
      pick[num_pick++] = n1;
      pick[num_pick++] = n2;
    }

  for(n1=m+1; n1<g->num_nodes; n1++)            /* add other nodes */
  {
    t=0;
    while(t<m)                                   /* insert m edges */
    {
      do
        n2 = (int) pick[(int) floor(drand48()*num_pick)];
      while(n2==n1);              /* chose pair of different nodes */
      if(!gs_edge_exists(g, n1, n2))
      {
	gs_insert_edge(g, n1, n2);
	pick[num_pick++] = n1;
	pick[num_pick++] = n2;
	t++;
      }
    }
  }
  free(pick);
}

/******************* graph_ausgabe() *********************/
/** Erstellt Ausgabedatei mit entsprechender            **/
/** Formatierung.                                       **/
/**                                                     **/
/** PARAMETERS: (*)= return-paramter                    **/
/**             (*)g: Graph                             **/
/**        num_nodes: Anzahl der Knoten                 **/
/**                                                     **/
/** RETURNS:                                            **/
/**     (nothing)                                       **/
/*********************************************************/
void graph_ausgabe(gs_graph_t *g, int num_nodes) {
  FILE *file;
  int i;
  file = fopen("graph_visu.dot", "w");
  fprintf(file, "graph scale_free_graph {\n");

  for(i=0; i<num_nodes; i++)
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

/******************* exportGraphDot() *********************/
/** Erstellt Ausgabedatei mit entsprechender             **/
/** Formatierung für die Visualisierung des Graphen      **/
/**                                                      **/
/** PARAMETERS: (*)= return-parameter                    **/
/**             (*)g: graph                              **/
/**                n: Knoten                             **/
/**                m: Anzahl der neuen Kanten            **/
/**                                                      **/
/**********************************************************/

void exportGraphDot(gs_graph_t *g,int n, int m)
{
    int i =0;
    FILE *file;
    char filename[1000];
    sprintf(filename, "Output/graph_visu_N%d_M%d.dot",n,m);
    file = fopen(filename, "w");
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

//TODO: print to console
void printDistances(double **dist, int n, int m)
{
    int i,j =0;
    FILE *file;
    char filename[1000];
    sprintf(filename, "Output_Normed/distance_N%d_M%d.dat",n,m);
    file = fopen(filename, "w");

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            fprintf(file,"%d, ", (int)dist[i][j]);
        }
        fprintf(file, "\n");
    }
}

void printHistogram(double *histogram, int n, int m)
{
    int i =0;
    FILE *file;
    char filename[1000];
    sprintf(filename, "Output_Normed/histogram_N%d_M%d.dat",n,m);
    file = fopen(filename, "w");

    for(i=1; i<n; i++)
    {
        fprintf(file,"%d %.10f\n", i, histogram[i]/(n*(n-1)/2));

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
                //TODO: inifintiy
                //dist[i][j] = g->num_nodes; maximaler Weg wenn edges in alle
                //Richtungen benutzt werden dürfen und Wegstrecke = 1
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

/**************runExperiments() **************/
/**                                         **/
/**                                         **/
/*********************************************/
void runExperiments(int runs, int n, int m)
{
    int i,j,k;
    gs_graph_t *g;
    double *histogram;
    histogram = (double*)malloc(sizeof(double)*n);

    double **dist;
    dist = (double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
    {
        dist[i] = (double*)malloc(sizeof(double)*n);
    }

    for (k=0;k<runs;k++)
    {
        g = gs_create_graph(n);
        gs_preferential_attachment(g, m);

        computeShortestPaths(g,dist);

        for (i=0;i<n;i++)
        {
            for (j=i+1;j<n;j++)
            {
                histogram[(int)dist[i][j]]++;
            }
        }
    }
    //TODO: free memory
    for (i=0;i<n;i++)
    {
        histogram[i] = histogram[i]/runs;
    }

    printHistogram(histogram,n,m);

}
