
#include "Snap.h"

int main(int argc, char* argv[])
{
	PUNGraph Graph = TUNGraph::New();

	Graph->AddNode(1);
	Graph->AddNode(3);
	Graph->AddNode(4);
	Graph->AddNode(5);
	Graph->AddEdge(1,5);
	Graph->AddEdge(1,3);
	Graph->AddEdge(1,4);
	TSnap::SaveEdgeList(Graph, "test.txt", "Save as tab-separated list of edges");

}