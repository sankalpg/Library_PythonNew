
#include "Snap.h"

static float  min_dist=1000000000, max_dist=-100000000000;

float ThshldArray[50] = {

0,
80,
180,
320,
500,
720,
980,
1280,
1620,
2000,
2420,
2880,
3380,
3920,
4500,
5120,
5780,
6480,
7220,
8000,
8820,
9680,
10580,
11520,
12500,
13520,
14580,
15680,
16820,
18000,
19220,
20480,
21780,
23120,
24500,
25920,
27380,
28880,
30420,
32000,
33620,
35280,
36980,
38720,
40500,
42320,
44180,
46080,
48020,
50001  
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
    
   
   createPatternGraphPerCollection(listFile, patternKNNExt, Graph, 0, ThshldArray[thresholdBin-1]);

    clusterCoff = TSnap::GetClustCf(Graph);
    
    sprintf(outFileName, "%d_ClusteringCoff.txt", thresholdBin);
    fp2 = fopen(outFileName, "w");
    fprintf(fp2, "%f\n", clusterCoff);
    fclose(fp2);
    
    sprintf(outInfoFile, "%d_NetworkInfo.txt", thresholdBin);
    TSnap::PrintInfo(Graph, "",outInfoFile, 0);

 
    //printf("Clustering coff is %f\n",clusterCoff);
    //printf("minimum and maximum distances are %f and %f respectively", min_dist, max_dist);

    /*TFOut FOut(outputNetworkFile); 
    Graph->Save(FOut);*/
    
    
    
    //TSnap::SaveEdgeList(Graph, outputNetworkFile, "Save as tab-separated list of edges");
    //TSnap::SavePajek(Graph, outputNetworkFile); 
	
	//TSnap::SaveEdgeList(Graph, "test.txt", "Save as tab-separated list of edges");
    
    //printf("Total number of files scanned = %d and correct ones are %d\n", ii, jj);

}
