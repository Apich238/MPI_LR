#include <mpi.h>
#include <malloc.h>
#include <stddef.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//using namespace::std;
int rank;
int size;

void ChannelLatency(int repeats) {
	double tt = 0;
	MPI_Status st;
	char* a = (char*)malloc(1);
	for (int t = 0; t < repeats; t++) {
		if (rank == 0) {
			for (int i = 1; i < size; i++) {
				double t0 = MPI_Wtime();
				MPI_Send(a, 0, MPI_BYTE, i, 1, MPI_COMM_WORLD);
				MPI_Recv(a, 0, MPI_BYTE, i, 2, MPI_COMM_WORLD, &st);
				tt += MPI_Wtime() - t0;
			}
		}
		else {
			MPI_Recv(a, 0, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &st);
			MPI_Send(a, 0, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
		}
	}
	tt /= 2.0* (size - 1)*repeats;
	if (rank == 0)
		printf("Channel latency is %2.8f seconds, found in %d repeats\n", tt, repeats);
}

void ChannelCapacity(int reps, int len) {
	double tt = 0;
	unsigned char* arr = (unsigned char*)malloc(len * sizeof(unsigned char));
	MPI_Status st;
	for (int j = 0; j < reps; j++) {
		if (rank == 0) {
			for (int i = 1; i < size; i++) {
				double a = MPI_Wtime();
				MPI_Send(arr, len, MPI_BYTE, i, 3, MPI_COMM_WORLD);
				MPI_Recv(arr, len, MPI_BYTE, i, 4, MPI_COMM_WORLD, &st);
				tt += MPI_Wtime() - a;
			}
		}
		else {
			MPI_Recv(arr, len, MPI_BYTE, 0, 3, MPI_COMM_WORLD, &st);
			MPI_Send(arr, len, MPI_BYTE, 0, 4, MPI_COMM_WORLD);
		}
	}
	free(arr);
	if (rank == 0) {
		double cap = 2.0*reps*len*(size - 1) / tt;
		printf("Channel capacity = %.0f bytes per second, found in %d repeats, data length = %d bytes\n", cap, reps, len);
	}
}

void TestBcast(int tsize, int datalength, int repeats) {
	MPI_Comm tmpcomm;
	MPI_Group all, tgroup;
	MPI_Comm_group(MPI_COMM_WORLD, &all);
	int* ranks = (int*)malloc(sizeof(int)*tsize);
	for (int i = 0; i < tsize; i++) ranks[i] = i;
	MPI_Group_incl(all, tsize, ranks, &tgroup);
	free(ranks);
	MPI_Comm_create(MPI_COMM_WORLD, tgroup, &tmpcomm);
	double bc_t = 0, bc_d, sr_t = 0, sr_d;
	if (tmpcomm != MPI_COMM_NULL) {
		char* buf = (char*)malloc(datalength * sizeof(char));
		int trank = 0;
		MPI_Comm_rank(tmpcomm, &trank);
		for (int q = 0; q < repeats; q++) {
			MPI_Barrier(tmpcomm);
			bc_d = MPI_Wtime();
			MPI_Bcast(buf, datalength, MPI_BYTE, 0, tmpcomm);
			bc_d = MPI_Wtime() - bc_d;
			bc_t += bc_d;
			MPI_Barrier(tmpcomm);
			sr_d = MPI_Wtime();
			if (trank == 0) {
				for (int i = 1; i < tsize; i++)
					MPI_Send(buf, datalength, MPI_BYTE, i, 5, tmpcomm);
			}
			else {
				MPI_Status st;
				MPI_Recv(buf, datalength, MPI_BYTE, 0, 5, tmpcomm, &st);
			}
			sr_d = MPI_Wtime() - sr_d;
			sr_t += sr_d;
		}
		MPI_Comm_free(&tmpcomm);
		free(buf);
	}
	MPI_Group_free(&tgroup);
	MPI_Group_free(&all);
	if (rank == 0)
		printf("Send&Recv time / Bcast time =%6.3f, found with data size %d bytes, excanges between %3d nodes in %5d repeats\n",
			sr_t / bc_t, datalength, tsize, repeats);
}

MPI_Status st;

#ifndef max
int max(int a,int b){return a>b?a:b;}
#endif

void MyReduce(char* my, char* reduced, int datalen, int redroot, MPI_Comm comm) {
	int myrank = 0, tsize = 0;
	MPI_Comm_rank(comm, &myrank);
	MPI_Comm_size(comm, &tsize);
	char* tmp1 = (char*)malloc(datalen * sizeof(char));
	char* tmp2 = (char*)malloc(datalen * sizeof(char));
	for (int i = 0; i < datalen; i++) tmp1[i] = my[i];
	int myrel = (myrank - redroot + tsize) % tsize;
	int pmax = 1;
	do {
		pmax *= 2;
	} while (pmax < tsize);
	for (int p = 1; p <= pmax; p *= 2) {
		int prel = myrel^p;
		bool isright = myrel > prel;
		if (max(myrel, prel) >= tsize) continue;
		int pabs = (prel + redroot) % tsize;
		if (isright) {
			MPI_Send(tmp1, datalen, MPI_BYTE, pabs, 6, comm);
			break;
		}
		else {
			MPI_Recv(tmp2, datalen, MPI_BYTE, pabs, 6, comm, &st);
			for (int i = 0; i < datalen; i++)
				tmp1[i] = tmp1[i] ^ tmp2[i];
		}
	}
	if (myrank == redroot)
		for (int i = 0; i < datalen; i++)
			reduced[i] = tmp1[i];
	free(tmp1);
	free(tmp2);
}

void TestReduce(int tsize, int datalength, int repeats) {
	MPI_Comm tmpcomm;
	MPI_Group all, tgroup;
	MPI_Comm_group(MPI_COMM_WORLD, &all);
	int* ranks = (int*)malloc(sizeof(int)*tsize);
	for (int i = 0; i < tsize; i++) ranks[i] = i;
	MPI_Group_incl(all, tsize, ranks, &tgroup);
	free(ranks);
	MPI_Comm_create(MPI_COMM_WORLD, tgroup, &tmpcomm);
	double r_t = 0, r_d, sr_t = 0, sr_d;
	if (tmpcomm != MPI_COMM_NULL) {
		char* my = (char*)malloc(datalength * sizeof(char)),
			*redused = (char*)malloc(datalength * sizeof(char));
		int trank = 0;
		MPI_Comm_rank(tmpcomm, &trank);
		int redroot = 1;
		char	*recieved = (char*)malloc(datalength * sizeof(char));
		MPI_Status st;
		for (int q = 0; q < repeats; q++) {
			MPI_Barrier(tmpcomm);
			r_d = MPI_Wtime();
			//usual reduce
			MPI_Reduce(my, redused, datalength, MPI_BYTE, MPI_BXOR, redroot, tmpcomm);
			r_d = MPI_Wtime() - r_d;
			r_t += r_d;
			MPI_Barrier(tmpcomm);
			sr_d = MPI_Wtime();
			//my reduce
			MyReduce(my, redused, datalength, redroot, tmpcomm);
			sr_d = MPI_Wtime() - sr_d;
			sr_t += sr_d;
		}
		free(recieved);
		MPI_Comm_free(&tmpcomm);
		free(my);
		free(redused);
	}
	MPI_Group_free(&tgroup);
	MPI_Group_free(&all);
	if (rank == 0)
		printf("Send&Recv time / Reduce time =%6.3f, found with data size %d bytes, excanges between %3d nodes in %5d repeats\n",
			sr_t / r_t, datalength, tsize, repeats);
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char* pend;
	double t = 0;
	/*int LatencyReps = argc > 1 ? strtol(argv[1], NULL, 10) : 1000000;
	int CapacityReps = argc > 2 ? strtol(argv[2], NULL, 10) : 100;
	int CapacityDataSize = argc > 3 ? strtol(argv[3], NULL, 10) : 256 * 1024 * 1024;
	int BCastReps = argc > 4 ? strtol(argv[4], NULL, 10) : 10;
	int BCastDataSize = argc > 5 ? strtol(argv[5], NULL, 10) : 256 * 1024 * 1024;
	int ReduceReps = argc > 6 ? strtol(argv[6], NULL, 10) : 10;
	int ReduceDataSize = argc > 7 ? strtol(argv[7], NULL, 10) : 512 * 1024 * 1024;*/
	int LatencyReps = 10000;
	////if (rank == 0) printf("Checking channels latency...\n");
	////t = MPI_Wtime();
	////ChannelLatency(LatencyReps);
	////t = MPI_Wtime() - t;
	////if (rank == 0) printf( "Latency measuring done in %f seconds\n", t);
	//=============================================================================================
	////int CapacityReps = 50;
	////int CapacityDataSize = 8 * 1024 * 1024;
	////if (rank == 0) printf("Checking channels capacity...\n");
	////for (int i = 1; i <=8; i++) {
	////	t = MPI_Wtime();
	////	ChannelCapacity(CapacityReps, CapacityDataSize*i);
	////	t = MPI_Wtime() - t;
	////}
	////if (rank == 0) printf("Capacity measuring done in %f seconds\n", t);
	//=============================================================================================
	////int BCastReps = 50;
	////int BCastDataSize = 8 * 1024 * 1024;
	////if (rank == 0) printf("Comparing Bcast and Send-Recv functions...\n");
	////for (int i = 1; i <=8; i++) {
	////	t = MPI_Wtime();
	////	for (int lim = 2; lim <= size; lim *= 2)
	////	{
	////		TestBcast(lim, BCastDataSize*i, BCastReps);
	////	}
	////	t = MPI_Wtime() - t;
	////	if (rank == 0) printf("Comparation done in %f seconds\n", t);
	////}
	//=============================================================================================
	int ReduceReps = 5;// 0;
	int ReduceDataSize = 8 * 1024 * 1024;
	if (rank == 0) printf("Comparing Reduce and Send-Recv functions...\n");
	for (int i = 1; i <= 8; i++) {
		t = MPI_Wtime();
		for (int lim = 2; lim <= size; lim *= 2)
		{
			TestReduce(lim, ReduceDataSize*i, ReduceReps);
		}
		t = MPI_Wtime() - t;
		if (rank == 0) printf("Comparation done in %f seconds\n", t);
	}
	//=============================================================================================
	if (rank == 0) printf("done\n");
	MPI_Finalize();
	return 0;
}