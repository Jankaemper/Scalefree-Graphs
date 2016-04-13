#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "graphs_lists.h"


/******************* gs_create_planar_graph() ********************/
/** Creates a graph with a fixed number of nodes 	    **/
/** in the [0,1]x[0,1] plane and with edges according   **/
/** to probability p_ij = f*[(1 + sqrt(N*pi)*d_ij) / α]^(−α) **/
/** PARAMETERS: (*)= return-paramter                    **/
/**      num_nodes: number of nodes                     **/
/**      f: constant for probability                    **/
/**      alpha: constant for probability                **/
/** RETURNS:                                            **/
/**     pointer to created graph                        **/
/*********************************************************/
gs_graph_t *gs_create_planar_graph(int num_nodes, double f, double alpha)
{
    srand48(time(NULL));
    gs_graph_t *g;
    int n,i,j;

    g = (gs_graph_t *) malloc(sizeof(gs_graph_t));  /* allocate */
    g->node = (gs_node_t *) malloc(num_nodes*sizeof(gs_node_t));

    g->num_nodes = num_nodes;                     	/* initialise nodes */
    for(n=0; n<num_nodes; n++)
    {
        g->node[n].neighbors = NULL;
        g->node[n].x_Coord = drand48();
        g->node[n].y_Coord = drand48();
    }

    double d_ij;
    for(i=0; i<num_nodes; i++)						/* init edges with defined probability */
    {
        for(j=i+1; j<num_nodes; j++)
        {
            d_ij=sqrt(pow((g->node[i].x_Coord - g->node[j].x_Coord),2)+pow((g->node[i].y_Coord - g->node[j].y_Coord),2));
            if (drand48() < f*pow(1/(1+sqrt(num_nodes*M_PI)*d_ij/alpha),alpha))
            {
                gs_insert_edge(g, i, j);
            }
        }
    }

    return(g);
}


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
    int max_pick;              /* maximum number of entries */

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

    for(n1=m+1; n1<g->num_nodes; n1++)   /* add other nodes */
    {
        t=0;
        while(t<m) /* insert m edges */
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

/******************* constant_probability() **************/
double *constant_propability(double constant,int current_nodes, int num_pick, int *pick)
{
    double *counter, *probability, proof;
    proof = 0.0;
    counter = (double *) malloc(num_pick*sizeof(double));
    probability = (double *) malloc(num_pick*sizeof(double));
    int i,j;

    for(j=0;j<current_nodes;j++)
    {
        for(i=0; i<num_pick; i++)
        {
            if (pick[i] == j)
            {
                counter[j] += 1.0;
            }

        }
        proof += (counter[j]+constant)/(num_pick+(constant*current_nodes));
        probability[j] = proof;
    }
    return(probability);
}



/******************* gs_preferential_attachment_constant() ********/
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
void gs_preferential_attachment_constant(gs_graph_t *g, int m, double constant)
{
    int i,t;
    int n1, n2;
    int *pick;            /* array which holds for each edge {n1,n2} */
    /* the numbers n1 and n2. Used for picking */
    /* nodes proportional to its current degree */
    int num_pick;              /* number of entries in 'pick' so far */
    int max_pick;                       /* maximum number of entries */
    double *probability, checkProbability;
    checkProbability = 0.0;
    if(g->num_nodes < m+1)
    {
        printf("graph too small to have at least %d edges per node!\n", m);
        exit(1);
    }
    max_pick = 2*m*g->num_nodes- m*(m+1);
    pick = (int *) malloc(max_pick*sizeof(int));
    probability = (double *)malloc(max_pick*sizeof(double));
    num_pick=0;

    for(n1=0; n1<m+1; n1++) /* start: complete subgraph of m+1 nodes */

        for(n2=n1+1; n2<m+1; n2++)
        {
            gs_insert_edge(g, n1, n2);
            pick[num_pick++] = n1;
            pick[num_pick++] = n2;
        }

    /*
       The nodes are added the same way like they are added in
       gs_preferential_attachment the only difference is the implementation.
       Instead of only using the 'pick' array to insert edgescreation,
       we calculate every probability for edges to be added to a node.
       This gives the possiblity to add a constant value to each possibility
       and therefore we are ableto influence the exponential structure of the edge
       distribution.
       */
    for(n1=m+1; n1<g->num_nodes; n1++)            /* add other nodes */
    {
        t=0;
        while(t<m)                                   /* insert m edges */
        {
            probability = constant_propability(constant, n1, num_pick, pick);
            double drand48();
            checkProbability = drand48();
            for(i=0;i<n1;i++)
            {
                //printf("%d %f %f\n",i, probability[i], checkProbability);
                if(checkProbability <= probability[i])
                {
                    n2 = i;
                    break;
                }
            }
            //printf("%d %d\n", n1, n2);
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
        while(next_neighb != NULL)
        {
            fprintf(file,"%d -- %d\n", i, next_neighb->info);
            next_neighb = next_neighb->next;
        }
    }
    fprintf(file, "}");
}

/**************printDistances() **************/
/** Writes the adjaceny matrix (distances) into a given file  **/
/** PARAMETERS: (*)= return-parameter                    **/
/**                dist: adajcency matrix                **/
/**                n: Nodes                              **/
/**                destination[]: Filename               **/
/**********************************************************/
void printDistances(double **dist, int n,char destination[])
{
    int i,j =0;
    FILE *file;
    file = fopen(destination, "w");

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            fprintf(file,"%f, ", dist[i][j]);
        }
        fprintf(file, "\n");
    }
}

/**************printedges() **************/
/** Writes the adjacency structure into a given File     **/
/** PARAMETERS: (*)= return-parameter                    **/
/**                g: given graph                        **/
/**                num_nodes: Nodes                      **/
/**                destination[]: Filename               **/
/*********************************************/
void printedges(gs_graph_t *g, int num_nodes, char destination[]) {
    FILE *file;
    int i;
    file = fopen(destination, "w");

    for(i=0; i<num_nodes; i++)
    {
        fprintf(file, "%f/%f: ",g->node[i].x_Coord,g->node[i].y_Coord);
        elem_t *next_neighb;
        next_neighb = g->node[i].neighbors;
        while(next_neighb != NULL)
        {
            fprintf(file,"%d: %f ", next_neighb->info, sqrt(pow(g->node[i].x_Coord - g->node[next_neighb->info].x_Coord,2)+pow(g->node[i].y_Coord - g->node[next_neighb->info].y_Coord,2)));
            next_neighb = next_neighb->next;
        }
        fprintf(file, "\n");
    }
}

/**************printHistogram() **************/
/** Writes the given histogram to the given destination file, each bin in a row    **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                histogram: given histogram              					       **/
/**                n: number of Nodes               					     	   **/
/**                destination[]: Filename        							       **/
/**                numBins: Histogram number of bins						       **/
/**                startIndex: first bin to be printed 						       **/
/*********************************************/
void printHistogram(double *histogram, int n, char destination[],int numBins, int startIndex)
{
    int i;
    FILE *file;
    file = fopen(destination, "w");

    for(i=startIndex; i<numBins; i++)
    {
        fprintf(file,"%d %.10f\n", i, histogram[i]);

    }
}

/**************printHistogramNormed() **************/
/** Writes the given histogram to the given destination file, each bin in a row    **/
/** Bin Values will be normed to overall probablity 1 							   **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                histogram: given histogram              					       **/
/**                n: number of Nodes               					     	   **/
/**                destination[]: Filename        							       **/
/**                numBins: Histogram number of bins						       **/
/**                startIndex: first bin to be printed 						       **/
/*********************************************/
void printHistogramNormed(double **histogram, int n, char destination[], int numBins, int startIndex, int runs)
{
    int i,j;
    FILE *file;
    file = fopen(destination, "w");
    double *mean, var, sum_mean;
    //counter number of histogram entries for the normation
    double columnSum = 0;
    double totalEntries =0;
    mean = (double*)malloc(sizeof(double)*numBins);

    for(i=startIndex; i<numBins; i++)
    {
        for (j = 0; j<runs;j++)
        {
            totalEntries += histogram[j][i];
        }
    }
    totalEntries /= runs;
            printf("%.10f\n", totalEntries);

    for(i=startIndex; i<numBins; i++)
    {
        columnSum = 0;
        for (j = 0; j<runs;j++)
        {
            histogram[j][i] = histogram[j][i]/totalEntries;
            columnSum += histogram[j][i];
        }
        mean[i] = columnSum/runs;
    }

    for(i=startIndex; i<numBins; i++)
    {
        var = 0;
        for (j = 0; j<runs;j++)
        {
            var += pow((mean[i]- histogram[j][i]),2);
        }
        var = sqrt(var/(runs-1));

        //if histogram labels start at 0, then add offset of 1 in order to avoid fitting problems (division by zero)
        if (startIndex == 0)
        {
            fprintf(file,"%d %.10f %.10f \n", i+1, mean[i], var);
        }
        else
        {
            fprintf(file,"%d %.10f %.10f \n", i, mean[i], var);
        }
    }
}

/**************computeShortestPaths() **************/
/** executes the Floyd Warshall Algorithm on the given Graph   				**/
/** PARAMETERS: (*)= return-parameter                 						**/
/**                g: given graph              					       		**/
/**                weighted: flag signaling if edges are weighted or not    **/
/*********************************************/
void gs_all_pair_shortest_paths(gs_graph_t *g,double **dist,int weighted)
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
                //workaround inifintiy
                dist[i][j] = 100000 ;
            }
        }
    }
    //set path distances for edges (if not weighted than 1 for all edges, otherwise euclidean distance)
    for(i=0; i<g->num_nodes; i++)
    {
        elem_t *next_neighb;
        next_neighb = g->node[i].neighbors;
        while(next_neighb != NULL)
        {
            //if weighted is false, then all edges should have weight 1
            if (weighted==0)
            {
                dist[i][next_neighb->info] = 1.0;
                dist[next_neighb->info][i] = 1.0;
            }
            //if weighted is true, then all edges should have their euclidean distance as length
            else
            {
                dist[i][next_neighb->info] = sqrt(pow(g->node[i].x_Coord - g->node[next_neighb->info].x_Coord,2)+pow(g->node[i].y_Coord - g->node[next_neighb->info].y_Coord,2));
                dist[next_neighb->info][i] = sqrt(pow(g->node[i].x_Coord - g->node[next_neighb->info].x_Coord,2)+pow(g->node[i].y_Coord - g->node[next_neighb->info].y_Coord,2));
            }
            next_neighb = next_neighb->next;
        }
    }
    //perform floid warshall algo
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

/**************fillHistogramDiscrete() **************/
/** for a given dist matrix with integer values   **/
/** packs the distances (representing shortest paths) into bins of histogram 	   **/
/** PARAMETERS: (*)= return-parameter          			   **/
/**                dist: distance matrix 		       **/
/**                histogram: given histogram  		       **/
/**                n: number of Nodes          	     	   **/
/*********************************************/
void fillHistogramDiscrete(double **dist, double ***histogram, int k, int n)
{
    int i,j;
    for (i=0;i<n;i++)
    {
        for (j=i+1;j<n;j++)
        {
            (*histogram)[k][(int)dist[i][j]]++;
        }
    }
}


/**************fillHistogramContinuous() **************/
/** for a given dist matrix and defined number of bins   **/
/** packs the distances (representing shortest paths) into bin of histogram. 	   **/
/** Non existing paths between vertices will be ignored							   **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                dist: distance matrix 			 						       **/
/**                histogram: given histogram              					       **/
/**                n: number of Nodes               					     	   **/
/**                numBins: number of Bins             					     	   **/
/*********************************************/
void fillHistogramContinuous(double **dist, double **histogram,int n, int numBins)
{
    int i,j;

    double minDist = n*sqrt(2)+1;
    double maxDist = 0;
    //find overall max and min distances
    for (i=0;i<n;i++)
    {
        for (j=i+1;j<n;j++)
        {
            if (dist[i][j] < minDist)
            {
                minDist = dist[i][j];
            }
            if (dist[i][j] > maxDist && dist[i][j] != 100000)
            {
                maxDist = dist[i][j];
            }
        }
    }

    //calc bin width
    double binWidth = (maxDist - minDist) / numBins;

    //assign all pair-wise distances from matrix to histogram
    for (i=0;i<n;i++)
    {
        for (j=i+1;j<n;j++)
        {
            if (dist[i][j] != 100000)
            {
                if (dist[i][j] == maxDist)
                {
                    (*histogram)[numBins-1]++;
                }
                else
                {
                    (*histogram)[(int)(floor((dist[i][j] - minDist) / binWidth))]++;
                }
            }
        }
    }
}

/**************computeMeanShortestPath() **************/
/** for a given dist matrix computes the mean value of   **/
/** all shortest paths in graph. Non existing connections between vertices will be ignored  **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                dist: distance matrix						       **/
/**                n: number of Nodes               					     	   **/
/*********************************************/
double computeMeanShortestPath(double **dist, int n)
{
    double sum=0;
    int i,j, countExisting =0;
    for (i=0;i<n;i++)
    {
        for (j=i+1;j<n;j++)
        {
            if (dist[i][j] != 100000)
            {
                countExisting ++;
                sum += dist[i][j];
            }
        }
    }
    return sum/countExisting;
}

/**************runExperimentsPlanar() **************/
/** creates scale-free graphs of fixed size and computes the histogram of  		   **/
/** all shortes paths in graph 					   						 		   **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                run: number of runs 			 						       	   **/
/**                n: number of Nodes               					     	   **/
/**                m: pref attachment parameter        					     	   **/
/*********************************************/
void runExperiments(int runs, int n, int m)
{
    int i,j,k;
    gs_graph_t *g;

    //init dist array
    double **dist;
    dist = (double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
    {
        dist[i] = (double*)malloc(sizeof(double)*n);
    }

    //perform runs
    double **histogram;
    histogram = (double**)malloc(sizeof(double*)*runs);
    for (i=0;i<runs;i++)
    {
        histogram[i] = (double*)malloc(sizeof(double)*n);
    }
    for (k=0;k<runs;k++)
    {
        //init scale free graph of size n
        g = gs_create_graph(n);
        gs_preferential_attachment(g, m);

        //compute all shortest paths with floyd warshall and
        gs_all_pair_shortest_paths(g,dist,0);

        //create histogram from distance matrix
        fillHistogramDiscrete(dist,&histogram,k,n);
    }

    //output in file
    char filename[1000];
    sprintf(filename, "OutputScaleFree_Normed/histogram_N%d_M%d.dat", n, m);
    printHistogramNormed(histogram, n, filename, n, 1, runs);

    free(histogram);
}


/**************runExperimentsPlanarConstant() **************/
/** creates scale-free graphs of fixed size and computes the histogram of  		   **/
/** all shortes paths in graph 					   						 		   **/
/** with constant value k_0 in probability **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                run: number of runs 			 						       	   **/
/**                n: number of Nodes               					     	   **/
/**                m: pref attachment parameter        					     	   **/
/*********************************************/
void runExperimentsConstant(int runs, int n, int m, double k_0)
{
    int i,j,k;
    gs_graph_t *g;


    //init dist array
    double **dist;
    dist = (double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
    {
        dist[i] = (double*)malloc(sizeof(double)*n);
    }

    //perform runs
    double *histogram;
    histogram = (double*)malloc(sizeof(double)*n);
    for (k=0;k<runs;k++)
    {
        //init scale free graph of size n
        g = gs_create_graph(n);
        gs_preferential_attachment_constant(g, m, k_0);

        //compute all shortest paths with floyd warshall and
        gs_all_pair_shortest_paths(g,dist,0);

        //create histogram from distance matrix
        //       fillHistogramDiscrete(dist,&histogram,n);
    }

    //output in file
    char filename[1000];
    sprintf(filename, "Output_constant/histogram_N%d_M%d_k%f.dat", n, m, k_0);
    //   printHistogramNormed(histogram, n, filename, n, 1);

    free(histogram);
}

/**************runExperimentsPlanar() **************/
/** creates graphs of fixed size on square [0,1]² and computes the histogram of  	**/
/** all shortes paths in graph 													 	**/
/** PARAMETERS: (*)= return-parameter                 							   	**/
/**                run: number of runs 			 						       	   	**/
/**                n: number of Nodes               					     	   	**/
/*********************************************/
void runExperimentsPlanar(int runs, int n)
{
    int i,j,k;
    gs_graph_t *g;

    //init dist array
    double **dist;
    dist = (double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
    {
        dist[i] = (double*)malloc(sizeof(double)*n);
    }

    //reasonable number of bins for histogram according to literature
    int numBins = sqrt(n);
    double *histogram = (double*)malloc(sizeof(double)*numBins);
    for (k=0;k<runs;k++)
    {
        //init planar graph evenly distributed in square [0,1]² of size n
        g = gs_create_planar_graph(n,1,3);

        //compute all shortest paths with floyd warshall
        gs_all_pair_shortest_paths(g,dist,1);

        //create histogram from distance matrix
        //     fillHistogramContinuous(dist,&histogram,n,numBins);
    }
    printDistances(dist,  n, "OutputPlanar_Normed/dists.dat");
    printedges(g, n,  "OutputPlanar_Normed/edges.dat");

    //output in file
    char filename[1000];
    sprintf(filename, "OutputPlanar_Normed/histogram_N%d.dat", n);
    //   printHistogramNormed(histogram, n, filename, numBins, 0);

    free(histogram);
}


/**************runMeanExperimentsPlanar() **************/
/** creates graph of varying sizes and writes              **/
/** mean length of shortest paths per graph size into file **/
/** PARAMETERS: (*)= return-parameter                 							   **/
/**                run: number of runs 			 						       	   **/
/**                stepSize: step size for iterating over graph sizes	     	   **/
/*********************************************/
void runMeanExperimentsPlanar(int runs, int stepSize)
{
    int i,j,k,n;
    gs_graph_t *g;
    FILE *file;
    file = fopen("OutputPlanar_Normed/meanValues.dat", "w");

    //write mean distance for different graph sizes into file
    for (n = 10; n <= 1000; n=n+stepSize)
    {
        //init array to 0
        double **dist = (double**)malloc(sizeof(double*)*n);
        for (i=0;i<n;i++)
        {
            dist[i] = (double*)malloc(sizeof(double)*n);
        }
        //perform runs
        double meanShortestPath= 0;
        for (k=0;k<runs;k++)
        {
            //init planar graph evenly distributed in square [0,1]² of size n
            g = gs_create_planar_graph(n,1,3);

            //compute all shortest paths with floyd warshall and write result in dist matrix
            gs_all_pair_shortest_paths(g,dist,1);

            //calc mean over all distances
            meanShortestPath += computeMeanShortestPath(dist,n);
        }

        meanShortestPath = meanShortestPath/runs;

        fprintf(file,"%d %f\n", n, meanShortestPath);
        printf("%d %f\n", n, meanShortestPath);
    }
}
