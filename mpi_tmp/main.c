#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

static const int MASTER_NODE_RANK = 0;

static const int INITIAL_MATRIX_TAG = 1;

static const int CURRENT_WORKER_TAG = 2;

int minimumEquationCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode);

int getEquationCountForDouble(int worldSize);

void fillMatrix(int numberOfEquations, int worldSize, int index, double *matrix);

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status status;
	clock_t start, end;
	double timeForNSendings, timeForNReduce;

	if (worldSize < 2) {
		printf("World is too small\n");
		return 1;
	}

	int isMaster = rank == MASTER_NODE_RANK;
	int equationCount;
	int equationsPerNode;
	int oneBlockSize;

	if (isMaster) {
		printf("World size is %d\n", worldSize);

		equationCount = getEquationCountForDouble(worldSize);
		equationsPerNode = equationCount / (worldSize - 1);
		oneBlockSize = (equationCount + 1) * equationsPerNode;
		printf("Equation count is %d\nwith %d per node\nwith size of %d\n", equationCount, equationsPerNode,
			oneBlockSize);
	}

	MPI_Bcast(&equationCount, 1, MPI_INT, MASTER_NODE_RANK, MPI_COMM_WORLD);
	MPI_Bcast(&equationsPerNode, 1, MPI_INT, MASTER_NODE_RANK, MPI_COMM_WORLD);
	MPI_Bcast(&oneBlockSize, 1, MPI_INT, MASTER_NODE_RANK, MPI_COMM_WORLD);

	double *matrix = malloc(oneBlockSize * sizeof(double));
	if (isMaster) {
		printf("Starting the sending of data\n");
		srand((unsigned int)time(0));

		int i;
		for (i = 1; i < worldSize; ++i) {
			fillMatrix(equationCount, worldSize, i - 1, matrix);
			MPI_Send(matrix, oneBlockSize, MPI_DOUBLE, i, INITIAL_MATRIX_TAG, MPI_COMM_WORLD);
		}
		free(matrix);
	}
	else {
		MPI_Recv(matrix, oneBlockSize, MPI_DOUBLE, MASTER_NODE_RANK, INITIAL_MATRIX_TAG, MPI_COMM_WORLD, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Данные разосланы.




	//    if (isMaster) {
	//        start = clock();
	//    }
	//
	//    char receive = 0;
	//    int sending;
	//    for (sending = 0; sending < numberOfSendings; ++sending) {
	//        MPI_Reduce(&rank, &receive, 1, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
	//    }
	//
	//    if (isMaster) {
	//        end = clock();
	//        timeForNReduce = (double) (end - start) / CLOCKS_PER_SEC;
	//        printf("time for %d reduce operations = %f \n", numberOfSendings, timeForNReduce);
	//        printf("OPS = %f\n", (double) numberOfSendings / timeForNReduce);
	//
	//    }
	if (!isMaster) {
		free(matrix);
	}
	MPI_Finalize();
}

int getEquationCountForDouble(int worldSize) {
	int doubleSize = sizeof(double);
	int oneHundreedMegaBytes = 100 * 1024 * 1024;
	int minValue = minimumEquationCount(doubleSize, worldSize - 1, oneHundreedMegaBytes);
	return minValue - minValue % (worldSize - 1) + worldSize - 1;
}

void fillMatrix(int numberOfEquations, int worldSize, int index, double *block) {
	int blockSize = numberOfEquations / (worldSize - 1);
	int i;
	int j;
	int mainIndex;
	for (i = 0; i < blockSize; ++i) {
		for (j = 0; j < numberOfEquations - 1; ++j) {
			mainIndex = index * (worldSize - 1) + i;
			block[i * numberOfEquations + j] = ((double)rand() / (double)RAND_MAX) * (j != mainIndex ? 10. : 100.);
		}
		block[i * numberOfEquations + numberOfEquations - 1] = ((double)rand() / (double)RAND_MAX) * 20 - 10;
	}
}

int minimumEquationCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
	return (int)ceil((sqrt((4. * minimumSizeOfSystemPerNode * numberOfComputationalNodes + dataTypeSizeInBytes)
		/ (double)dataTypeSizeInBytes) - 1.)
		/ 2.);
}