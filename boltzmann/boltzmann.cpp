#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

using namespace std;


int rank, size;

double frand(double min, double max) {
	return min + rand()*((max - min) / RAND_MAX);
}

double** InitGrid(int l) {
	double** buf = (double**)malloc(sizeof(double*)*l);
	for(int i=0;i<l;i++)
		buf[i] = (double*)malloc(sizeof(double)*l);
	srand(0);
	for (int i = 0; i < l; i++)
		for (int j = 0; j < l; j++)
			buf[i][j] = frand(0, 5);
	return buf;
}

void FreeGrid(double** g,int l) {
	for (int i = 0; i < l; i++) free(g[i]);
	free(g);
}

//           N
//     NW    /\    NE
//           ||
//  W <======()======> E
//           ||
//     SW    \/    SE
//           S
//
// Y
// /\
// ||
// ()=>X


struct CalcPtsGroup2d {
	CalcPtsGroup2d *N,*S,*W,*E,*NW,*NE,*SW,*SE;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
};


int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int gridlen = 8;
	double** grid = InitGrid(gridlen);


	


	FreeGrid(grid,gridlen);
	MPI_Finalize();
}

