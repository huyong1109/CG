/*******************************************************************
 *       Filename:  cross.c                                     
 *                                                                 
 *    Description:        N(=9) processes compute Jacobi method of                                 
 *                          A (9*9) for M(=1) times.                                 
 * 			Simplified from a application. 
 *        Version:  1.0                                            
 *        Created:  03/31/2012 08:34:45 PM                                 
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
#include"mpi.h"
#define SUB_SIZE 3
#define N_PROCESS 9
#define CHANGEDIM(AA,xx,yy)  *((double*)(AA)+(SUB_SIZE+2)*(xx)+(yy))

/* Cross
 * passing massage to up, down, left ,right 
 * */
int Cross(double **A, int myid,int subsize,int np,MPI_Status *status ){
	int left, right,up, down;
	if (myid /3 == 0)
		up = MPI_PROC_NULL;
	else 
		up = myid -3;
	if (myid /3 == 2)
		down = MPI_PROC_NULL;
	else 
		down = myid +3;
	if (myid %3 == 0)
		left = MPI_PROC_NULL;
	else 
		left = myid -1;
	if (myid %3 == 2)
		right = MPI_PROC_NULL;
	else 
		right = myid +1;
/* send massage to up proc and receive from down proc */
	int i,j;
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid ==0 )
		printf("Send massage to up and Recv from down\n");	

	MPI_Sendrecv(&CHANGEDIM(A,1,1),subsize,MPI_DOUBLE,up,1,&CHANGEDIM(A,subsize+1,1),subsize,MPI_DOUBLE,down,1,MPI_COMM_WORLD,status);
		if(up != MPI_PROC_NULL){
		printf(" proc %d send [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,1,i+1));}
		printf("] to    %d\n",up);}
		if(down != MPI_PROC_NULL ){
		printf(" proc %d get  [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,subsize+1,i+1));}
		printf("] from  %d\n",down);}

/* send massage to down proc and receive from up proc */	
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid ==0 )
		printf("Send massage to down and Recv from up\n");	

	MPI_Sendrecv(&CHANGEDIM(A,subsize,1),subsize,MPI_DOUBLE,down,2,&CHANGEDIM(A,0,1),subsize,MPI_DOUBLE,up,2,MPI_COMM_WORLD,status);
		if(down != MPI_PROC_NULL){
		printf(" proc %d send [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,subsize,i+1));}
		printf("] to    %d\n",down);}
		if(up != MPI_PROC_NULL ){
		printf(" proc %d get  [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,0,i+1));}
		printf("] from  %d\n",up);}

/* send massage to left proc and receive from right  proc */	
	MPI_Datatype A_COL;
	MPI_Type_vector(subsize,1,subsize+2,MPI_DOUBLE,&A_COL);
	MPI_Type_commit(&A_COL);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid ==0 )
		printf("Send massage to left and Recv from right\n");	

	MPI_Sendrecv(&CHANGEDIM(A,1,1),1,A_COL,left,3,&CHANGEDIM(A,1,subsize+1),1,A_COL,right,3,MPI_COMM_WORLD,status);
		if(left != MPI_PROC_NULL){
		printf(" proc %d send [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,i+1,1));}
		printf("] to    %d\n",left);}
		if(right != MPI_PROC_NULL ){
		printf(" proc %d get  [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,i+1,subsize+1));}
		printf("] from  %d\n",right);}
/* send massage to down proc and receive from up proc */	
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid ==0 )
		printf("Send massage to left and Recv from right\n");	
	
	MPI_Sendrecv(&CHANGEDIM(A,1,subsize),1,A_COL,right,4,&CHANGEDIM(A,1,0),1,A_COL,left,4,MPI_COMM_WORLD,status);
		if(right != MPI_PROC_NULL){
		printf(" proc %d send [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,i+1,subsize));}
		printf("] to    %d\n",right);}
		if(left != MPI_PROC_NULL ){
		printf(" proc %d get  [",myid);
	for (i =0;i<subsize; i++){
		printf("%f\t",CHANGEDIM(A,i+1,0));}
		printf("] from  %d\n",left);}
	
	return 0;
}

/* Jacobi processing 
 * 	A[i][j]= 0.25(A[i-1][j]+A[i+1][j]+A[i][j+1]+A[i][j-1])
 * */
int Jacobi(double **A, int myid,int subsize,int np){
	int i,j;
	double B[subsize][subsize];
	for(i=1;i<subsize+1;i++)
		for(j=1;j<subsize+1;j++)
		{
			B[i-1][j-1]=0.25*(CHANGEDIM(A,i-1,j)+ CHANGEDIM(A,i+1,j)+CHANGEDIM(A,i,j+1)+CHANGEDIM(A,i,j-1));
		}
	for(i=1;i<subsize+1;i++){
		for(j=1;j<subsize+1;j++){
		CHANGEDIM(A,i,j)=B[i-1][j-1];}}
	return 0;
}


/* CrossShow
 * print the value of each sub matrix, with extra side numbers 
 * */
int CrossShow(double **A, int myid,int subsize)
{
	int i,j,k;
	printf("Value in proc %d\n",myid);
		for(j=0;j<subsize+2;j++){
			for(k=0;k<subsize+2;k++)
			printf("%1.3f\t",CHANGEDIM(A,j,k));
			printf("\n\n");}
	
	return 0;
}
/* main function
 * */
int main(int argc, char **argv)
{
	int subsize = SUB_SIZE;
	double A[subsize+2][subsize+2];
	int count, tag, myid,np, i,j,k;
	MPI_Status status;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	/* check processes  */
	if(np < N_PROCESS){
		printf("!!!%d processes is needed!!!\n",N_PROCESS );
		return -1;
	}	
	/*  set initiate ,all get 0.0 */
	for(i=0;i<subsize+2;i++)
		for(j=0;j<subsize+2;j++)
			A[i][j]= 0.0;

	/* set four sides element to 8.0 */ 
	for(i =0;i<3;i++){
		if(myid == i){
		for(j=1;j<subsize+1;j++)
		A[1][j]= 8.0;
		}
	}
	for(i =6;i<9;i++){
		if(myid == i){
		for(j=1;j<subsize+1;j++)
		A[subsize][j]= 8.0;
		}	
	}
	for(i =0;i<9;i=i+3){
		if(myid == i){
		for(j=1;j<subsize+1;j++)
		A[j][1]= 8.0;}
	}
	for(i =2;i<9;i=i+3){
		if(myid == i){
		for(j=1;j<subsize+1;j++)
		A[j][subsize]= 8.0;}
	}
	if (myid ==0){	
	printf("----------Initiate Matrix  --------\n");}
	for(i =0;i<9;i++){
		if(myid ==i ){
		CrossShow(A,myid,subsize);}
		MPI_Barrier(MPI_COMM_WORLD);	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid ==0){	
	printf("----------Passing massage -----begin---\n");}	
	Cross(A,myid, subsize,np,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
	if (myid ==0){	
	printf("----------Passing massage -----end---\n");	
	printf("----------Matrix show after  Passage --------\n");}
	for(i =0;i<9;i++){
		if(myid ==i ){
		CrossShow(A,myid,subsize);}
		MPI_Barrier(MPI_COMM_WORLD);	
	}
	Jacobi(A, myid, subsize,np);
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid ==0){	
	printf("----------Matrix show after  Jacobi --------\n");}
	for(i =0;i<9;i++){
		if(myid ==i ){
		CrossShow(A,myid,subsize);}
		MPI_Barrier(MPI_COMM_WORLD);	
	}
	MPI_Finalize();
	return 0;
}
