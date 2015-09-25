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
  gs_graph_t *g;
  int num_real, m;                                /* number of realizations */
  int argz = 1;                   /* for treating command line arguments */

  if(argc != 3)
  {
    printf("USAGE %s <N> <m> \n", argv[0]);
    exit(1);
  }

  num_real = 1000;
  num_nodes = atoi(argv[argz++]);
  m = atoi(argv[argz++])               /* read number of nodes */
  //  sscanf(argv[argz++], "%lf", &c2);     /* read inter-group connectivity */

  g = gs_create_graph(num_nodes);
  gs_preferential_attachment(g, m)



  return(0);

}
