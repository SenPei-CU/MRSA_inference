  //Agent-based modeling of MRSA transmission.
/*********************************************************************
 * agentbasedsimulation.c
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 * Adapted from the code by Shawn Lankton (http://www.shawnlankton.com/2008/03/getting-started-with-mex-a-short-tutorial/)
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlmxhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *mxnl, *mxpart, *mxstate, *mxpara, *mxnodes, *mxappearnodes, *mxcol_intro, *mxcol_envir, *mxalpha, *mxpk, *mxmu, *mxwardcnt, *mxwardcap, *mxwardcontam, *mxnewstate, *mxstat, *mxnewinf, *mxnewwardcontam, *mxnewcol_envir, *mxisinnet;
    //mxstat: new colonized in hospitals; new transmitted infection from colonized in hospitals;
    //new infection from colonized outside hospitals
    const mwSize *dims;
    double *nl, *part, *state, *para, *nodes, *appearnodes, *col_intro, *col_envir, *alpha, *pk, *mu, *wardcnt, *wardcap, *wardcontam, *newstate, *stat, *newinf, *newwardcontam, *newcol_envir, *isinnet;
    int dimx, dimy, Nmax, num_para, num_nodes, num_app, nllength, num_ward;

//associate inputs
    mxnl = mxDuplicateArray(prhs[0]);
    mxpart = mxDuplicateArray(prhs[1]);
    mxstate = mxDuplicateArray(prhs[2]);
    mxpara = mxDuplicateArray(prhs[3]);
    mxnodes = mxDuplicateArray(prhs[4]);
    mxappearnodes = mxDuplicateArray(prhs[5]);
    mxcol_intro = mxDuplicateArray(prhs[6]);
    mxcol_envir = mxDuplicateArray(prhs[7]);
    mxalpha = mxDuplicateArray(prhs[8]);
    mxpk = mxDuplicateArray(prhs[9]);
    mxmu = mxDuplicateArray(prhs[10]);
    mxwardcnt = mxDuplicateArray(prhs[11]);
    mxwardcap = mxDuplicateArray(prhs[12]);
    mxwardcontam = mxDuplicateArray(prhs[13]);

//figure out dimensions
    dims = mxGetDimensions(prhs[2]);//state
    dimy = (int)dims[0]; dimx = (int)dims[1]; Nmax = dimy;
    dims = mxGetDimensions(prhs[3]);//para:beta,I0
    num_para = (int)dims[0];
    dims = mxGetDimensions(prhs[4]);//nodes
    num_nodes = (int)dims[0];
    dims = mxGetDimensions(prhs[5]);//appearnodes
    num_app = (int)dims[0];
    dims = mxGetDimensions(prhs[0]);//nl
    nllength = (int)dims[0];
    dims = mxGetDimensions(prhs[12]);//wardcap
    num_ward = (int)dims[0];

//associate outputs
    mxnewstate = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    mxstat = plhs[1] = mxCreateDoubleMatrix(5,1,mxREAL);
    mxnewinf = plhs[2] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    mxnewwardcontam = plhs[3] = mxCreateDoubleMatrix(num_ward,1,mxREAL);
    mxnewcol_envir = plhs[4] = mxCreateDoubleMatrix(Nmax,1,mxREAL);
    
//create variables
    mxisinnet = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    nl = mxGetPr(mxnl);
    part = mxGetPr(mxpart);
    state = mxGetPr(mxstate);
    nodes = mxGetPr(mxnodes);
    appearnodes = mxGetPr(mxappearnodes);
    col_intro = mxGetPr(mxcol_intro);
    col_envir = mxGetPr(mxcol_envir);
    alpha = mxGetPr(mxalpha);
    pk = mxGetPr(mxpk);
    mu = mxGetPr(mxmu);
    wardcnt = mxGetPr(mxwardcnt);
    wardcap = mxGetPr(mxwardcap);
    wardcontam = mxGetPr(mxwardcontam);
    para = mxGetPr(mxpara);//beta,I0
    newstate = mxGetPr(mxnewstate);
    newinf = mxGetPr(mxnewinf);
    stat = mxGetPr(mxstat);
    stat[0] = 0; //new colonized due to transmission
    stat[1] = 0; //new colonized due to environmental contamination
    stat[2] = 0; //new infection from colonized outside hospitals
    stat[3] = 0; //new infection from transmittd colonization
    stat[4] = 0; //new infection from environmental colonization
    newwardcontam = mxGetPr(mxnewwardcontam);
    newcol_envir = mxGetPr(mxnewcol_envir);
    
    isinnet = mxGetPr(mxisinnet);
//do something
    int i, j, node, nei, t, T = 7;
    double v;
    double beta, I0, theta, lambda, alphai, pi, mui;
    //parameters
    beta = para[0]; I0 = para[1]; theta = para[3]; lambda = para[4];
//     printf("%f,%f,%f,%f\n",beta,I0,theta,lambda);
    for (i = 0; i < num_nodes; i++)
        isinnet[(int)(nodes[i]-1)] = 1;
    //////////////////////////////////////
    ////////////////////////
    //starting T looping
    for (t = 1; t <= T; t++){
    //copy state to newstate
    for (i = 0; i < Nmax; i++)
        newstate[i] = state[i];
    //copy wardcontam to newwardcontam, decay
    for (i = 0; i < num_ward; i++)
        newwardcontam[i] = (1-lambda)*wardcontam[i];
    //copy col_envir to newcol_envir
    for (i = 0; i < Nmax; i++)
        newcol_envir[i] = col_envir[i];
    //update nodes in the network
    //state: 0 - susceptible, 1 - colonized, 2 - infected    
    for (i = 0; i < num_nodes; i++) {
        node = (int)(nodes[i]-1);
        //get the room number
        int ward = wardcnt[i+(t-1)*num_nodes]; 
        if (state[node] == 0) {//susceptible  
//             //calculate roommate number other than node
//             double ni=0;
//             for (j = part[node]-1; j < part[node+1]-1; j++){
//                 if (nl[j+t*nllength] > 0)
//                     ni++;
//             }
            //neighbors
            for (j = part[node]-1; j < part[node+1]-1; j++){
                nei = (int)(nl[j]-1);
                if (nl[j+t*nllength] > 0){//connected at t
                    if (state[nei] > 0) {//colonized or infected
                        v = (double)rand()/RAND_MAX;
                        int ni = wardcap[ward-1];//get room capacity
                        if (ni < 2)
                            ni = 2;
//                         printf("%f\n",beta/(ni-1));
                        if (v < beta/(ni-1))
                            newstate[node] = 1;
                    }
                }
            }
            if (newstate[node] == 1)
                stat[0]++;
            if (ward>0){//currently in a ward
                //environmental contaminination
                double eta = wardcontam[ward-1];//contamination rate
                v = (double)rand()/RAND_MAX;
                if (v < eta){
                    newstate[node] = 1;
                    stat[1]++;
                    newcol_envir[node] = 1;
                }
            }
        }
        if (state[node] == 1) {//colonized
            alphai = alpha[node];
            pi = pk[node]*alphai;
            v = (double)rand()/RAND_MAX;
            if (v < alphai)//decolonzied
                newstate[node] = 0;
            else if (v > 1-pi)//infection
                newstate[node] = 2;
            if (newstate[node] == 2){//infected
                newinf[node] = 1;
                if (col_envir[node]>0)//colonization due to environmental
                    stat[4]++;
                else if (col_intro[node]>0)//colonization outside hospitals
                    stat[2]++;
                else
                    stat[3]++;//colonization due to transmission
            }
            //contribute to environmental contamination
            if ((newstate[node]!=0)&&(ward>0)){//currently in a ward
                int ni = wardcap[ward-1];//get room capacity
                if (ni < 1)
                    ni = 1;
                newwardcontam[ward-1] += theta/ni;
            }
        }
        if (state[node] == 2) {//infection
            mui = mu[node];
            v = (double)rand()/RAND_MAX;
            if (v < mui)//recoverd
                newstate[node] = 0;
            //contribute to environmental contamination
            if ((newstate[node]!=0)&&(ward>0)){//currently in a ward
                int ni = wardcap[ward-1];//get room capacity
                if (ni < 1)
                    ni = 1;
                newwardcontam[ward-1] += theta/ni;
            }
        }
    }
  //update nodes ever appear that outside hospitals
    for (i = 0; i < num_app; i++) {
        node = (int)(appearnodes[i]-1);
        if (isinnet[node] == 0) {//not in the network
            if (state[node] == 1) {//colonized
                alphai = alpha[node];
                pi = pk[node]*alphai;
                v = (double)rand()/RAND_MAX;
                if (v < alphai)//decolonzied
                    newstate[node] = 0;
                else if (v > 1-pi)//infection
                    newstate[node] = 2;
            }
            if (state[node] == 2) {//infection
                mui = mu[node];
                v = (double)rand()/RAND_MAX;
                if (v < mui)//recoverd
                    newstate[node] = 0;
            }
        }
    }
    //copy newstate to state
    for (i = 0; i < Nmax; i++)
        state[i] = newstate[i];
    //copy newwardcontam to wardcontam
    for (i = 0; i < num_ward; i++)
        wardcontam[i] = newwardcontam[i];
    //copy newcol_envir to col_envir
    for (i = 0; i < Nmax; i++)
        col_envir[i] = newcol_envir[i];
    }//end t loop
    return;
}
