
#include "Snap.h"

static float  min_dist=1000000000, max_dist=-100000000000;
  
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


int computeClusterCoffCurve(char *listFile, char *patternKNNExt, , char *ClusterCoffFile)
{

}

int main(int argc, char* argv[])
{
    FILE *fp1, *fp2;
    char *listFile, *patternKNNExt,  tempFilename[400]= {'\0'}, *outputNetworkFile, patternKNNFile[400]= {'\0'};
    int verbos;
    int ii=0, jj=0;
    long long int ind1, ind2;
    float dist;
    double clusterCoff;
    
    PUNGraph Graph = TUNGraph::New();
    
    
    if(argc < 4 || argc > 5)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    
    patternKNNExt = argv[1];
    listFile = argv[2];
    outputNetworkFile = argv[3];     
    if( argc == 5 ){verbos = atoi(argv[4]);}
    
   
   createPatternGraphPerCollection(listFile, patternKNNExt, Graph, 0,10000);

    clusterCoff = TSnap::GetClustCf(Graph);
    
    printf("Clustering coff is %f\n",clusterCoff);
    printf("minimum and maximum distances are %f and %f respectively", min_dist, max_dist);

    /*TFOut FOut(outputNetworkFile); 
    Graph->Save(FOut);*/
    
    
    
    //TSnap::SaveEdgeList(Graph, outputNetworkFile, "Save as tab-separated list of edges");
    TSnap::SavePajek(Graph, outputNetworkFile); 
	
	//TSnap::SaveEdgeList(Graph, "test.txt", "Save as tab-separated list of edges");
    
    //printf("Total number of files scanned = %d and correct ones are %d\n", ii, jj);

}
