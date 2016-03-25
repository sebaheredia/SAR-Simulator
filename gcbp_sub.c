/*
 * =============================================================
 * gcbp_sub.c 
 * Matlab mex file
 * im = gcbp_sub(im,cbuf,ptnew,g_off,dx,dy,wt,fo,ro,c,pi,rmin,rmax,dr,nc,nx,ny)
 * Convolution Backprojection sub routine for use in gcbp.m
 *
 * =============================================================
 */

/* $Revision: 1.0 $ */
#include <math.h>
#include "mex.h"

/* Computational subroutine */
void gcbp_sub(double *imr,double *imi,double *cbufr, double *cbufi,double *ptnew,
	    double *g_off,double dx, double dy,double wt,double fo,double ro,
 	    double c,double pi,double rmin, double rmax,double dr,int nc,int lf,
 	    int nx,int ny,int ncb)
{
  int i,j,k,count=0;
  double vt[3];
  double vtmag,rr,alpha,xxx,grid3,grid2,grid1;
  double t,p,rs,q,f;
  int n1,n2;
  double cvalr = 0.0;
  double cvali = 0.0;
  double cvalr2,cvali2;
  
  /* Perform the loops of the complex vectors. */
  grid3 = 0.0;
  for (i = 0; i < ny; i++) {
      grid2 = (ny/2 - i )*dy + *(g_off+1);
      for (j = 0; j < nx; j++) {
        grid1 = (j - nx/2 )*dx + *(g_off+0);
      
        vt[0] = *(ptnew) - grid1;
        vt[1] = *(ptnew+1) - grid2;
        vt[2] = *(ptnew+2) - grid3;
		
	    vtmag = vt[0]*vt[0] + vt[1]*vt[1] + vt[2]*vt[2];
	    vtmag = sqrt(vtmag);	
	    
	    rr = ro;	
	    alpha = 2.0*(vtmag-rr);
	    xxx = alpha/2.0;
	    
	    if (xxx <= rmax && xxx >= rmin){
		        t = alpha/c;
            	p = 2*pi*fo*t;
            	rs = nc/2+(alpha/2)/dr;   /*distance along pulse in samples*/
            	n1 = rs;
            	n2 = n1+1;
            	q = rs-n1;   /*fraction away from n1*/

/*if samples lie within cbuf array, then backproject (interpolate)*/
if (lf == 0){
    /*use linear interpolation*/
            if (n1 >= 1 && n2 <= nc){
            
                cvalr = ((1.0-q)**(cbufr+n1) + q**(cbufr+n2))*cos(p) - 
                        ((1.0-q)**(cbufi+n1) + q**(cbufi+n2))*sin(p);
                
                cvali = ((1.0-q)**(cbufr+n1) + q**(cbufr+n2))*sin(p) +
                        ((1.0-q)**(cbufi+n1) + q**(cbufi+n2))*cos(p);
                
                *(imr+count) = *(imr+count) + cvalr*wt;
                *(imi+count) = *(imi+count) + cvali*wt; 
		    }
            else{
                *(imr+count) = 0.0;
                *(imi+count) = 0.0;
            }

        count++;
}
else{
            if ( n1 >= (lf+1) && n1 <= (nc-lf) ){
                    cvalr = 0.0;
                    cvali = 0.0;
                    for (k = -lf; k <= lf; k++){
                        f = sin(pi*(k-q) + 1.0e-10)/(pi*(k-q) + 1.0e-10);
                        cvalr = cvalr + f**(cbufr+n1+k);
                        cvali = cvali + f**(cbufi+n1+k);
                    }
                
                    cvalr2 = cvalr*cos(p) - cvali*sin(p);
                    cvali2 = cvalr*sin(p) + cvali*cos(p);
                
                    *(imr+count) = *(imr+count) + cvalr2*wt;
                    *(imi+count) = *(imi+count) + cvali2*wt; 
            }
            else{
                    *(imr+count) = 0.0;
                    *(imi+count) = 0.0;
            }
        count++;
}
    }
   }
}
}
/* The gateway routine. */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  double  *imr,*imi,*cbufr,*cbufi,*ptnew,*g_off,dx,dy,wt,fo,ro,c,pi,rmin,rmax,dr;
  int     nx,ny,ncb,nc,lf;

    
  /* Check for the proper number of arguments. */
  if (nrhs != 16)
    mexErrMsgTxt("16 inputs required.");
  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");

  /* Get pointers to real and imaginary parts of the input. */
  imr = mxGetPr(prhs[0]);
  imi = mxGetPi(prhs[0]);
  cbufr = mxGetPr(prhs[1]);
  cbufi = mxGetPi(prhs[1]);
  ptnew = mxGetPr(prhs[2]);
  g_off = mxGetPr(prhs[3]);
  
  /* Get the scalar input x. */
  dx = mxGetScalar(prhs[4]);
  dy = mxGetScalar(prhs[5]);
  wt = mxGetScalar(prhs[6]);
  fo = mxGetScalar(prhs[7]);
  ro = mxGetScalar(prhs[8]);
  c = mxGetScalar(prhs[9]);
  pi = mxGetScalar(prhs[10]);
  rmin = mxGetScalar(prhs[11]);
  rmax = mxGetScalar(prhs[12]);
  dr = mxGetScalar(prhs[13]);
  nc = mxGetScalar(prhs[14]);
  lf = mxGetScalar(prhs[15]);
  
  /* Get the dimensions of the matrix input im. */
  nx = mxGetM(prhs[0]);
  ny = mxGetN(prhs[0]);

  
  /* Get the length of input vector. */
  ncb = mxGetM(prhs[1]);

  /* Create a new array and set the output pointer to it. */
  plhs[0] = mxCreateDoubleMatrix(nx, ny, mxCOMPLEX);
  imr = mxGetPr(plhs[0]);
  imi = mxGetPi(plhs[0]);
  
  /* Call the C subroutine. */
  gcbp_sub(imr,imi,cbufr,cbufi,ptnew,g_off,dx,dy,wt,fo,ro,c,pi,rmin,rmax,dr,nc,lf,nx,ny,ncb);

  return;
}
