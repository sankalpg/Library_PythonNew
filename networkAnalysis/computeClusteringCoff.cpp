
#include "Snap.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


static float  min_dist=1000000000, max_dist=-100000000000;

float ThshldArray[50] = {
0.0,
 208.247,
 832.986,
1874.219,
3331.945,
5206.164,
7496.876,
10204.082,
13327.780,
16867.972,
20824.656,
25197.834,
29987.505,
35193.669,
40816.327,
46855.477,
53311.120,
60183.257,
67471.887,
75177.010,
83298.626,
91836.735,
100791.337,
110162.432,
119950.021,
130154.102,
140774.677,
151811.745,
163265.306,
175135.360,
187421.908,
200124.948,
213244.481,
226780.508,
240733.028,
255102.041,
269887.547,
285089.546,
300708.038,
316743.024,
333194.502,
350062.474,
367346.939,
385047.897,
403165.348,
421699.292,
440649.729,
460016.660,
479800.083,
500001.000
};


  
int createPatternGraphPerFile(char *fileName, PUNGraph Graph, double dist_min, double dist_max)
{
	FILE *fp2;
	long long int ind1, ind2;
	int ii=0;
	float dist;

	fp2 = fopen(fileName, "r");
        if (fp2==NULL)
        {
            printf("Error opening file %s\n", fileName);
        }
        else
        {
            while(fscanf(fp2, "%lld\t%lld\t%f\n",&ind1,&ind2, &dist)!=EOF)
            {
                if ((dist >= dist_max)||(dist <= dist_min))
                {
                    continue;    
                }
                if (!Graph->IsNode(ind1)) Graph->AddNode(ind1);
                if (!Graph->IsNode(ind2)) Graph->AddNode(ind2);
                if (!Graph->IsEdge(ind1, ind2)) Graph->AddEdge(ind1, ind2);

            }
            fclose(fp2);
        }
        
     return 1;
  
}

int createPatternGraphPerCollection(char *listFile, char *patternKNNExt, PUNGraph Graph, double dist_min, double dist_max)
{
    FILE *fp1;
    int ii=0;
    char  tempFilename[400]= {'\0'},  patternKNNFile[400]= {'\0'};

    fp1 = fopen(listFile, "r");
	if (fp1==NULL)
    {
        printf("Error opening file %s\n", listFile);
        return 0;
    }
        
    while(fgets(tempFilename, 400, fp1))
    {
        sscanf(tempFilename, "%[^\n]s\n", patternKNNFile);
        
        strcat(patternKNNFile,patternKNNExt);
        
        printf("processing file %d\t%s\n", ii+1,patternKNNFile);
	    createPatternGraphPerFile(patternKNNFile, Graph, dist_min, dist_max);
        memset(tempFilename, '\0', sizeof(char)*400);
        memset(patternKNNFile, '\0', sizeof(char)*400);
        ii++;
    }
    fclose(fp1);
    
    return 1;
 
}
/*
 * THis function randomizes a network. This was needed to show how a current network clusters better compared to a random network.
 * NOTE that during the process of randomization we preserve several network topological characteristics
 * we preserve, number of notes, number of edges, in and out degree of every node.
 */
int randomizeGraph(PUNGraph Graph, int NSwapsMultiple)
{
    long int nNodes = Graph->GetNodes();
    long int nSwaps = NSwapsMultiple*nNodes;
    long int cnt = nSwaps, nChanges=0;
    int nId1, nId2, nId1nn, nId2nn, nn1, nn2, rnn1, rnn2;
    srand (time(NULL));
    TUNGraph::TNodeI node1, node2, node1nn, node2nn;
    //printf("Hello1\n");
    while(cnt >0)
    {
        nChanges++;
        //printf("-----------------------CNT %d----------------------------\n",cnt);
        //get a random note
        nId1 = Graph->GetRndNId();
        nId2 = Graph->GetRndNId();
        //printf("nodes are %d and %d\n",nId1, nId2);
        //Check if there is a link between these two nodes, if yes, just skip the following steps
        if(Graph->IsEdge(nId1, nId2)){continue;}
        
        
        node1 = Graph->GetNI(nId1);
        nn1 = node1.GetOutDeg();
        //printf("Out degree1 %d\n",nn1);
        rnn1 = rand()%nn1;
        //printf("random number 1 = %d\n",rnn1);
        nId1nn = node1.GetNbrNId(rnn1);
        //printf("Hello3\n");
        if(Graph->IsEdge(nId1nn, nId2)){continue;}
        
        node2 = Graph->GetNI(nId2);
        nn2 = node2.GetOutDeg();
        //printf("Out degree2 %d\n",nn2);
        rnn2 = rand()%nn2;
        //printf("random number 2 = %d\n",rnn2);
        nId2nn = node2.GetNbrNId(rnn2);
        //printf("Hello4\n");
        if(Graph->IsEdge(nId1, nId2nn)){continue;}
        
        if(Graph->IsEdge(nId1nn, nId2nn)){continue;}
        
        //printf("Hello5\n");
        Graph->AddEdge(nId1, nId2);
        Graph->AddEdge(nId1nn, nId2nn);
        
        //printf("Hello6\n");
        
        Graph->DelEdge(nId1, nId1nn);
        Graph->DelEdge(nId2, nId2nn);
        //rintf("Hello7\n");
        
        cnt--;
    }
    //rintf("total iterations %d\n",nChanges);
    
    return 1;
}


/*int computeClusterCoffCurve(char *listFile, char *patternKNNExt, , char *ClusterCoffFile)
{

}*/

int main(int argc, char* argv[])
{
    FILE *fp1, *fp2;
    char *listFile, *patternKNNExt,  tempFilename[400]= {'\0'}, *outputNetworkFile, patternKNNFile[400]= {'\0'}, outFileName[400]= {'\0'}, outInfoFile[400]= {'\0'};;
    int verbos;
    int ii=0, jj=0;
    long long int ind1, ind2;
    float dist;
    double clusterCoff;
    int thresholdBin;
    
    PUNGraph Graph = TUNGraph::New();
    
    
    if(argc < 5 || argc > 6)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    patternKNNExt = argv[1];
    listFile = argv[2];
    outputNetworkFile = argv[3];
    thresholdBin = atoi(argv[4]);
    if( argc == 6 ){verbos = atoi(argv[5]);}
    
   
   createPatternGraphPerCollection(listFile, patternKNNExt, Graph, 0, 10.0*ThshldArray[thresholdBin-1]);

    clusterCoff = TSnap::GetClustCf(Graph);
    
    sprintf(outFileName, "%d_ClusteringCoff.txt", thresholdBin);
    fp2 = fopen(outFileName, "w");
    fprintf(fp2, "%f\n", clusterCoff);
    fclose(fp2);
    
    sprintf(outInfoFile, "%d_NetworkInfo.txt", thresholdBin);
    //TSnap::PrintInfo(Graph, "",outInfoFile, 0);
    
    
    //save graph
    TSnap::SaveEdgeList(Graph, "GraphEdgeListFile", "Save as tab-separated list of edges");
    TSnap::SavePajek(Graph, "GraphPajekFile"); 
    
    printf("Number of nodes are: %d\n",Graph->GetNodes());
    randomizeGraph(Graph, (int)1);
    
    clusterCoff = TSnap::GetClustCf(Graph);
    
    sprintf(outFileName, "%d_ClusteringCoff_RANDOM.txt", thresholdBin);
    fp2 = fopen(outFileName, "w");
    fprintf(fp2, "%f\n", clusterCoff);
    fclose(fp2);
    
    sprintf(outInfoFile, "%d_NetworkInfo_RANDOM.txt", thresholdBin);
    //TSnap::PrintInfo(Graph, "",outInfoFile, 0);
    
    
    

 
    //printf("Clustering coff is %f\n",clusterCoff);
    //printf("minimum and maximum distances are %f and %f respectively", min_dist, max_dist);

    /*TFOut FOut(outputNetworkFile); 
    Graph->Save(FOut);*/
    
    
    

	
	//TSnap::SaveEdgeList(Graph, "test.txt", "Save as tab-separated list of edges");
    
    //printf("Total number of files scanned = %d and correct ones are %d\n", ii, jj);

}
