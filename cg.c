/*******************************************************************
 *       Filename: cg.c                                     
 *                                                                 
 *    Description:  CG solve equations      A*X = b
 *        Version:  1.0                                            
 *        Created:  03/31/2012 08:34:45PM                                 
 *       Revision:  none                                           
 *       Compiler:  gcc                                           
 *                                                                 
 *         Author:  Hu Yong                                      
 *          Email:  huyong1109@gmail.com                                        
 *        Company:  HPC tsinghua                                      
 *                                                                 
 *******************************************************************/
#include<stdio.h>
#include<stdlib.h>
#define N 1000

int main(int argc, char *argv )
{

       double A[N][N];
       double x[N],b[N];
       readpro(A,x,b);
}
    /*
     * read data
     */
int readpro(double **A, double *x, double *b) 
{
       int DSIZE = sizeof(double);
       FILE *Abin;
       FILE *bbin;
       FILE *xbin;
       Abin = fopen("A.bin","rb");
       if (!Abin)
       {
           printf("File A.bin not exist!!!");
       }
       bbin = fopen("b.bin","rb");
       if (!bbin)
       {
           printf("File b.bin not exist!!!");
       }
       xbin = fopen("x.bin","rb");
       if (!xbin)
       {
           printf("File x.bin not exist!!!");
       }

           printf("Begin read !!!");
       fread(A, DSIZE, N*N,Abin); 
       fread(b, DSIZE, N,bbin); 
       fread(x, DSIZE, N,xbin); 
       fclose(Abin);
       fclose(xbin);
       fclose(bbin);
}
