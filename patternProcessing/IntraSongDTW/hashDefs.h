
#ifndef HASHDEFS_H

#define HASHDEFS_H

//typdefs
#define DATATYPE       double
#define DISTTYPE        double
#define INDTYPE        long long
#define PID_DEFAULT1    -1
#define PID_DEFAULT2    -2
#define PID_DEFAULT3    -3
#define FID_DEFAULT1    -1


// other Constants
#define INF FLT_MAX
#define LOG2  0.693147180559945
#define EPS 0.0000000000000001
#define MAXNTEMPOFACTORS 10



// inlined functions
//#define computeLBkimFL(a,b,c,d) ((a-b)*(a-b)) + ((c-d)*(c-d))
//#define oneDtwoD(i,j,N) (i*N + j)

//Extern declarations
extern const int combAllwd_5[5][5];
extern const int combAllwd_3[3][3];


#endif //HASHDEFS_H

