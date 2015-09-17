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


int main(int argc, char **argv)
{
  int num_nodes;                             /* number of nodes in graph */
  double c1, c2;                               /* average connectivities */
  gs_graph_t *g;
  int *size;                                          /* of each cluster */
  double largest;                                        /* largest size */
  int *group;             /* for each node, IDF of its preassigned group */
  int *id;                           /* for each node, id of its cluster */
  int num_cl;                                      /* number of clusters */
  int t, cl;                                            /* loop counters */
  int num_real;                                /* number of realizations */
  int argz = 1;                   /* for treating command line arguments */
  double size0, size1, size2,size4;  
                                /* moments of sizes of largest component */
  int num_percol;

  if(argc != 3)
  {
    printf("USAGE %s <N> <c2>\n", argv[0]);
    exit(1);
  }

  num_real = 1000;
  num_nodes = atoi(argv[argz++]);                /* read number of nodes */
  sscanf(argv[argz++], "%lf", &c2);     /* read inter-group connectivity */

  g = gs_create_graph(num_nodes);
  group = (int *) malloc(num_nodes*sizeof(int));
  id = (int *) malloc(num_nodes*sizeof(int));
  /*gs_insert_edge(g, 0, 5);
  gs_insert_edge(g, 5, 7);
  gs_insert_edge(g, 7, 0);
  gs_insert_edge(g, 3, 8);
  gs_write_graph(g, stdout);
  */
  for(t=0; t<num_nodes; t++)
    group[t] = (2*t)/num_nodes;

  
  for(c1=0.2; c1<2.0; c1+=0.1)        /* loop over different connectivities */
  {
    size0=0.0; size1=0.0; size2=0.0; size4=0.0;       /* reset statistics */
    num_percol = 0;
    for(t=0; t<num_real; t++)
    {
      gs_clear_graph(g);                        /* generate random graph */
      gs_block_model(g, 2*c1/num_nodes, 2*c2/num_nodes, group);
      size = gs_clusters(g, id, &num_cl);
      largest=size[0];                     /* look for largest component */
      for(cl=1; cl<num_cl; cl++)
	if(size[cl] > largest)
	  largest = size[cl];
      size0++;                                          /* do statistics */
      largest /= num_nodes;
      if(largest > 0.3)
	num_percol++;
      size1 += largest;
      largest = largest*largest;
      size2 += largest;
      largest = largest*largest;
      size4 += largest;
      free(size);
    }
    printf("%f %f %f %f %f  %f %d\n", c1, c2, size1/size0, 
	   sqrt( (size2/size0-size1*size1/(size0*size0))/num_real), 
	   0.5*(3-(size4)/(size2*size2/size0)), (double) num_percol/size0, 
	   num_nodes);

  }

 
  gs_delete_graph(g);

  free(group);
  free(id);
  return(0);

}
