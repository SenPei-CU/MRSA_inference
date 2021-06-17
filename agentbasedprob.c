//Agent-based modeling of MRSA transmission.
/*********************************************************************
 * agentbasedsimulation.cpp
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
    mxArray *mxnl, *mxpart, *mxS, *mxC, *mxpara, *mxnodes, *mxappearnodes, *mxalpha, *mxpk, *mxwardcnt, *mxwardcap, *mxwardcontam, *mxnewS, *mxnewC, *mxnewwardcontam, *mxisinnet;
    
    const mwSize *dims;
    double *nl, *part, *S, *C, *para, *nodes, *appearnodes, *alpha, *pk, *wardcnt, *wardcap, *wardcontam, *newS, *newC, *newwardcontam, *isinnet;
    int dimx, dimy, Nmax, num_para, num_nodes, num_app, nllength, num_ward;

//associate inputs
    mxnl = mxDuplicateArray(prhs[0]);
    mxpart = mxDuplicateArray(prhs[1]);
    mxS = mxDuplicateArray(prhs[2]);
    mxC = mxDuplicateArray(prhs[3]);
    mxpara = mxDuplicateArray(prhs[4]);
    mxnodes = mxDuplicateArray(prhs[5]);
    mxappearnodes = mxDuplicateArray(prhs[6]);
    mxalpha = mxDuplicateArray(prhs[7]);
    mxpk = mxDuplicateArray(prhs[8]);
    mxwardcnt = mxDuplicateArray(prhs[9]);
    mxwardcap = mxDuplicateArray(prhs[10]);
    mxwardcontam = mxDuplicateArray(prhs[11]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[2]);//state
    dimy = (int)dims[0]; dimx = (int)dims[1]; Nmax = dimy;
    dims = mxGetDimensions(prhs[4]);//para:beta,I0
    num_para = (int)dims[0];
    dims = mxGetDimensions(prhs[5]);//nodes
    num_nodes = (int)dims[0];
    dims = mxGetDimensions(prhs[6]);//appearnodes
    num_app = (int)dims[0];
    dims = mxGetDimensions(prhs[0]);//nl
    nllength = (int)dims[0];
    dims = mxGetDimensions(prhs[10]);//wardcap
    num_ward = (int)dims[0];

//associate outputs
    mxnewS = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    mxnewC = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    mxnewwardcontam = plhs[2] = mxCreateDoubleMatrix(num_ward,1,mxREAL);
//create variables
    mxisinnet = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    nl = mxGetPr(mxnl);
    part = mxGetPr(mxpart);
    S = mxGetPr(mxS);
    C = mxGetPr(mxC);
    nodes = mxGetPr(mxnodes);
    appearnodes = mxGetPr(mxappearnodes);
    alpha = mxGetPr(mxalpha);
    pk = mxGetPr(mxpk);
    wardcnt = mxGetPr(mxwardcnt);
    wardcap = mxGetPr(mxwardcap);
    wardcontam = mxGetPr(mxwardcontam);
    para = mxGetPr(mxpara);//beta,I0
    newS = mxGetPr(mxnewS);
    newC = mxGetPr(mxnewC);
    newwardcontam = mxGetPr(mxnewwardcontam);
    
    isinnet = mxGetPr(mxisinnet);
//do something
    int i, j, node, nei, t, T = 7;
    double v;
    double beta, I0, theta, lambda, alphai;
    beta = para[0]; I0 = para[1]; theta = para[3]; lambda = para[4];
    double decol, newcol, treat, inf;
    for (i = 0; i < Nmax; i++)
        isinnet[i] = 0;
    for (i = 0; i < num_nodes; i++){
        isinnet[(int)(nodes[i]-1)] = 1;
    }
    ////////////////////
    for (t = 1; t <= T; t++){
        //copy state to newstate
        for (i = 0; i < Nmax; i++){
            newS[i] = S[i];
            newC[i] = C[i];
        }
        //copy wardcontam to newwardcontam, decay
        for (i = 0; i < num_ward; i++){
            newwardcontam[i] = (1-lambda)*wardcontam[i];
        }
        //update nodes in the network

        //state: 0 - susceptible, 1 - colonized, 2 - infected    
        for (i = 0; i < num_nodes; i++) {
            node = (int)(nodes[i]-1);
            //get the room number
            int ward = wardcnt[i+(t-1)*num_nodes]; 
            //////////prepare parameters
            alphai = alpha[node];
            ///////////////
            decol = alphai*C[node];
            newcol = 0;
            for (j = part[node]-1; j < part[node+1]-1; j++){
                nei = (int)(nl[j]-1);
                if (nl[j+t*nllength] > 0){//connected at t
                    int ni = wardcap[ward-1];//get room capacity
                    if (ni < 2)
                        ni = 2;
                    newcol += beta/(ni-1)*S[node]*C[nei];
                }
            }
            /////////////add environmental contamination
            if (ward>0)//currently in a ward
                newcol += wardcontam[ward-1];
            //////////////
            newS[node] += decol-newcol;
            newC[node] += newcol-decol;
            ////////////update wardcontam
            if (ward>0){//currently in a ward
                int ni = wardcap[ward-1];//get room capacity
                if (ni < 2)
                    ni = 2;
                newwardcontam[ward-1] += theta/(ni)*newC[node];
            }
        }
        //update nodes ever appear that outside hospitals
        for (i = 0; i < num_app; i++) {
            node = (int)(appearnodes[i]-1);
            if (isinnet[node] == 0) {//not in the network
                //////////prepare parameters
                alphai = alpha[node];
                //////////
                decol = alphai*C[node];
                ////////
                newS[node] += decol;
                newC[node] += -decol;
            }
        }
        //copy newstate to state
        for (i = 0; i < Nmax; i++){
            S[i] = newS[i];
            C[i] = newC[i];
        }
        //copy newwardcontam to wardcontam
        for (i = 0; i < num_ward; i++)
            wardcontam[i] = newwardcontam[i];
    }
    return;
}
