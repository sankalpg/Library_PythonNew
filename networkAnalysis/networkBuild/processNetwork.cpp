
#include "Snap.h"

int main(int argc, char* argv[])
{
    FILE *fp1, *fp2;
    char *graphFile;    
    int verbos;
    int ii=0, jj=0;
    long long int ind1, ind2;
    float dist;
    
    if(argc < 2 || argc > 3)
    {
        printf("\nInvalid number of arguments!!!\n");
        exit(1);
    }
    graphFile = argv[1];
    if( argc == 3 ){verbos = atoi(argv[2]);}
    
    printf("Trying to read %s\n",graphFile);
    TFIn FIn(graphFile); 
    PNGraph Graph = TNGraph::Load(FIn);

}