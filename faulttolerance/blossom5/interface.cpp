#include <stdio.h>
#include "PerfectMatching.h"
#include "GEOM/GeomPerfectMatching.h"



extern "C" {

int* blossom_match(int node_num, int edge_num, int nodes1[], int nodes2[], int weights[])
{

	struct PerfectMatching::Options options;

	char* filename = NULL;
	bool check_perfect_matching = false;
	int* edges;

	PerfectMatching *pm = new PerfectMatching(node_num,edge_num);

	for(int k=0; k<edge_num;++k)
	{
		pm->AddEdge(nodes1[k],nodes2[k],weights[k]);
	}
	
	options.verbose = false;
	pm->options = options;
	pm->Solve();

	int i, j;
	int* partners = new int[node_num];
	for (i=0; i<node_num; i++)
	{
		j = pm->GetMatch(i);
		partners[i]=j;
	}

	delete pm;

	return partners;
}

}
