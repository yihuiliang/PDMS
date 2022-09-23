#include <math.h>
#include "mex.h"


/*
    paretomember returns the logical Pareto membership of a set of points.

    synopsis:  front = paretofront(objMat)

    created by Yi Cao
    
    y.cao@cranfield.ac.uk
    
    for compiling type 

    mex paretofront.c
    
*/


void FDMO(bool * front, double * M, unsigned int row, unsigned int col);

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	bool * front;
	double * M;
	unsigned int row, col;
	const int  *dims;
    
	if(nrhs == 0 || nlhs > 1)
	{
	    printf("\nsynopsis:   front = FDMO(X)");
	    plhs[0]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
	    return;
	}
	
	M = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	row = dims[0];
	col = dims[1];
	
	
	
	/* ----- output ----- */

	plhs[0]    = mxCreateLogicalMatrix (row , 1);
	front = (bool *) mxGetPr(plhs[0]);
	
	
	/* main call */
	FDMO(front,  M, row, col);
}

void FDMO(bool * front, double * M, unsigned int row, unsigned int col)
{
	
	//init
    unsigned int *candidatelist;    
    candidatelist = (unsigned int *)mxMalloc(row*sizeof(unsigned int));	
    for(int i=0;i<row;i++)
		candidatelist[i] = i;
	
	unsigned int *ptr_lb, *ptr_ub, *ptr_cmp, *ptr_base;
	ptr_base = candidatelist;
	ptr_lb = candidatelist;
	ptr_ub = candidatelist+row-1;
	ptr_cmp = ptr_lb+1;
	// cycle
	unsigned int tmp;
	unsigned int j,j1,j2;
	while(ptr_lb!=ptr_ub)
	{
		while(ptr_cmp-1!= ptr_ub)
		{
			//cmpare A and B
			bool flag1 = false,flag2 = false;
			for (j=0,j1=*ptr_lb,j2=*ptr_cmp;j<2;j++) {
				if (M[j1] > M[j2]) {
					flag1 = true;
				}else{
					flag2 = true;
				}
                if(j == 1)
                    break;
                j1+=row;
                j2+=row;
			}
			if(flag1)
			{
				if(flag2){
					//////////////////////A as good as B///////////////////
					ptr_cmp++;
				}
				else{
					/////////////////B better than A//////////////////////
					tmp = *ptr_lb;
					*ptr_lb = *ptr_cmp;
					*ptr_cmp = *ptr_ub;
					*ptr_ub = tmp;
					ptr_ub--;
					ptr_cmp = ptr_lb+1;
				}
			}
			else
			{
				/////////////////A better than B//////////////////////
				//swap ptr_ub and ptr_cmp
				tmp = *ptr_cmp;
				*ptr_cmp = *ptr_ub;
				*ptr_ub = tmp;
				ptr_ub--;
			}

		}
		if(ptr_lb==ptr_ub)
			break;
		ptr_lb ++;
		ptr_cmp = ptr_lb+1;				
	}		
	for(unsigned int *ptr = ptr_base;ptr<=ptr_lb;ptr++)
		front[*ptr] = true;
    mxFree(candidatelist); 
}