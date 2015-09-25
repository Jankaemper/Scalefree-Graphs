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

/*
Jedem neuen Knoten sollen m Kanten hinzugefügt werden. D.h. der erste Schritte
ist zu überprüfen ob die gewünschte Anzahl an Knoten überhaupt möglich ist.
Beispielsweise gilt, wenn der Ausgangsgraph aus 5 Knoten besteht, kann der neue
Knoten maximal 5 Kanten haben.
Um diese Bedingung zu erfüllen folgende Abfrage:
*/

  if(g->num_nodes < m+1)
  {
      printf("graph too small to have at least %d edges per node!\n", m);
      exit(1);
  }

  // --------------------------------------------------------------------------
  /*
  max_pick bestimmt die größe des Arrays: 2 * Anzahl der Kanten die hinzugefügt
  werden sollen * die Anzahl der Knoten - (Kanten*(Kanten+1)).
  Für einen Ursrpungsgraph von 5 Knoten und m=2 gilt also:
  2*2*5 - (2*(2+1)) = 20 - 6 = 14
  */

  max_pick = 2*m*g->num_nodes- m*(m+1);
  pick = (int *) malloc(max_pick*sizeof(int));
  num_pick=0;

/* Es werden soviele Kanten wie möglich erstellt. Kompletter subgraph mit
neuen Knoten. Pick array wird hochgezählt um die Kanten zu speichern.
*/

  for(n1=0; n1<m+1; n1++) /* start: complete subgraph of m+1 nodes */

    for(n2=n1+1; n2<m+1; n2++)
    {
      gs_insert_edge(g, n1, n2);
      pick[num_pick++] = n1;
      pick[num_pick++] = n2;
    }
// ----------------------------------------------------------------------------

/*
Schleife über alle alten Knoten. n2 wird zufällig gewählt und mit n1, welches
die alten Knoten abläuft verglichen. Der Wertebereich von n2 entspricht
der gespeicherten Werte im pick-array.
Wenn nun n2 = n1 ist wird geprüft ob bereits eine Kante vorhanden ist,
sonst wird sie hinzugefügt. Je mehr
*/

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
