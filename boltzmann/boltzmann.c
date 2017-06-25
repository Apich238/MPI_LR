#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <memory.h>

#define LATTICE_DIRECTIONS 9

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

const double elementalVectors[LATTICE_DIRECTIONS][2] = {{0,  0},
                                                        {1,  0},
                                                        {0,  1},
                                                        {-1, 0},
                                                        {0,  -1},
                                                        {1,  1},
                                                        {-1, 1},
                                                        {-1, -1},
                                                        {1,  -1}};

typedef struct {
    double density;//���������������� ���������
    double velocity[2];    //���������������� ��������, 0 - �������������, 1 - �����������
} MacroNode;

typedef struct {
    double particleDistribution[LATTICE_DIRECTIONS];        //������������� ������ �� ������������
    double tmp[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //������� �����
    double relaxationTime, latticeSpeed;        //����� ���������� � �������� �����
    GridNode **nodes;
} Grid;

typedef struct {
    int first;
    int last;
} RowBounds;

void testVectorFunctions();

void testVectorSum();

void testVectorMultiply();

void testVectorModulus();

void testScalarMultiplication();

void testModelFunctions();

void testDensity();

void testVelocity();

void boundsComputationTest();

void sumVector(double *first, double *second, double *result);

void multiplyVector(double *vector, double multiplier, double *result);

double modulusOfVector(double *vector);

double scalarMultiplication(double *first, double *second);

double cosBetweenVectors(double *first, double *second);

RowBounds getMyBounds(int gridWidth, int computationalProcessorsCount, int index);

double calculateDensity(double *particleDistribution);

void calculateVelocity(double *particleDistribution, double macroscopicDensity, double latticeSpeed, double *result);

#pragma region vectorF

void sumVector(double *first, double *second, double *result) {
	int i;
	for (i = 0; i < 2; ++i) {
		result[i] = first[i] + second[i];
	}
}

void multiplyVector(double *vector, double multiplier, double *result) {
	int i;
	for (i = 0; i < 2; ++i) {
		result[i] = vector[i] * multiplier;
	}
}

double modulusOfVector(double *vector) {
	return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

double scalarMultiplication(double *first, double *second) {
	double result = 0;
	int i;
	for (i = 0; i < 2; ++i) {
		result += first[i] * second[i];
	}
	return result;
}

double cosBetweenVectors(double *first, double *second) {
	return scalarMultiplication(first, second) / (modulusOfVector(first) * modulusOfVector(second));
}
#pragma endregion
#pragma region tests

void testModelFunctions() {
	testDensity();
	testVelocity();
	boundsComputationTest();
}

void boundsComputationTest() {
	RowBounds firstBounds = getMyBounds(5, 3, 0);
	if (firstBounds.first != 0 && firstBounds.last != 1) {
		exit(103);
	}
	RowBounds secondBounds = getMyBounds(5, 3, 1);
	if (secondBounds.first != 2 && secondBounds.last != 3) {
		exit(104);
	}
	RowBounds thirdBounds = getMyBounds(5, 3, 2);
	if (thirdBounds.first != 4 && thirdBounds.last != 4) {
		exit(105);
	}

	firstBounds = getMyBounds(4, 3, 0);
	if (firstBounds.first != 0 && firstBounds.last != 1) {
		exit(106);
	}
	secondBounds = getMyBounds(4, 3, 1);
	if (secondBounds.first != 2 && secondBounds.last != 2) {
		exit(107);
	}
	thirdBounds = getMyBounds(4, 3, 2);
	if (thirdBounds.first != 3 && thirdBounds.last != 3) {
		exit(108);
	}
}

void testVelocity() {
	double particleDistribution[LATTICE_DIRECTIONS] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	double density = calculateDensity(particleDistribution);
	double velocity[2];
	calculateVelocity(particleDistribution, density, 0.1, velocity);
	if (velocity[0] != -0.0055555555555555601 || velocity[1] != -0.016666666666666666) {
		exit(102);
	}
}

void testDensity() {
	double particleDistribution[LATTICE_DIRECTIONS] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	double density = calculateDensity(particleDistribution);
	if (density != (8 + 0) / 2 * 9) {
		exit(101);
	}
}

void testScalarMultiplication() {
	double first[2] = { 1, 2 };
	double second[2] = { 3, 4 };
	double scalar = scalarMultiplication(first, second);
	if (scalar != 1 * 3 + 2 * 4) {
		exit(5);
	}
}

void testVectorModulus() {
	double first[2] = { 4, 3 };
	double modulus = modulusOfVector(first);
	if (modulus != 5.) {
		exit(4);
	}
}

void testVectorMultiply() {
	double first[2] = { 1, 2 };
	double result[2];
	multiplyVector(first, 2, result);
	if (result[0] != 2 || result[1] != 4) {
		exit(3);
	}
}

void testVectorSum() {
	double first[2] = { 1, 2 };
	double second[2] = { 3, 4 };
	double result[2];
	sumVector(first, second, result);
	if (result[0] != 4 || result[1] != 6) {
		exit(2);
	}
}

void testVectorFunctions() {
	testVectorSum();
	testVectorMultiply();
	testVectorModulus();
	testScalarMultiplication();
}

#pragma endregion
#pragma region bounds
/**
* @param gridWidth ������ ���������� �����
* @param computationalProcessorsCount ���������� ������������
* @param index ����� ����������� ������� � 0
* @return ������� ������ � ��������� ����� � �����.
*/
RowBounds getMyBounds(int gridWidth, int computationalProcessorsCount, int index) {
	int remainder = gridWidth % computationalProcessorsCount;
	RowBounds res;
	res.first = gridWidth / computationalProcessorsCount * index + (index < remainder ? index : remainder);
	res.last = gridWidth / computationalProcessorsCount * (index + 1) - 1 + (index < remainder ? index + 1 : remainder);
	return res;
}

int minimumRowCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
	return (int)ceil((sqrt((minimumSizeOfSystemPerNode * numberOfComputationalNodes)
		/ dataTypeSizeInBytes)));
}
#pragma endregion
#pragma region LBM specific
void fillTempFieldForNode(const Grid *grid, const GridNode *upperBound, int hasUpperBound, const GridNode *lowerBound,
                          int hasLowerBound, int row, int column);
/**
 * @param particleDistribution ������������� ������ �� ������������
 * @param macroscopicDensity ���������������� ��������� � �����
 * @param latticeSpeed �������� �����
 * @param result ���������������� �������� � �����
 */
void calculateVelocity(double *particleDistribution, double macroscopicDensity, double latticeSpeed, double *result) {
    double temp[2];
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = 0;
    }
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        multiplyVector((double *) elementalVectors[direction], particleDistribution[direction], temp);
        multiplyVector((double *) temp, latticeSpeed, temp);
        sumVector(result, temp, result);
    }

    multiplyVector(result, 1. / macroscopicDensity, result);
}

/**
 * @param particleDistribution ������������� ������ �� ������������
 * @return ���������������� ��������� � �����
 */
double calculateDensity(double *particleDistribution) {
    double density = 0;
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        density += particleDistribution[direction];
    }
    return density;
}

/**
 * @param direction �����������
 * @param latticeSpeed �������� �����
 * @param velocity ���������������� ��������
 * @return ����������� ��� ���������� ������������ ������������� �� ������������
 */
double s(int direction, double latticeSpeed, double *velocity) {
    double scalar = scalarMultiplication((double *) elementalVectors[direction], velocity);
    return weights[direction] * 3 *
           (scalar + (3 * pow(scalar, 2) - scalarMultiplication(velocity, velocity)) / (latticeSpeed * 2)) /
           latticeSpeed;
}

/**
 * @param latticeSpeed �������� �����
 * @param density ���������������� ���������
 * @param velocity ���������������� ��������
 * @param result ����������� ������������� �� ������������ (OUT)
 */
void calculateEquilibriumDistribution(double latticeSpeed, double density, double *velocity, double *result) {
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = (1 + s(direction, latticeSpeed, velocity)) * density * weights[direction];
    }
}
void Streaming(Grid *pg, int rank, int worldSize) {
	int hasUpperBound = rank != 1;
	int hasLowerBound = rank != (worldSize - 1);
	GridNode *upperBound, *lowerBound;
	size_t rowSize = sizeof(GridNode) * pg->width;


	//����� �������� �������� �����
	if (hasUpperBound) {
		upperBound = malloc(rowSize);
		//�������� ��, ��� ����� ��������
		memcpy(upperBound, pg->nodes[0], rowSize);
	}
	if (hasLowerBound) {
		lowerBound = malloc(rowSize);
		//�������� ��, ��� ����� ��������
		memcpy(lowerBound, pg->nodes[pg->height - 1], rowSize);
	}
	MPI_Status status;
	for (int i = 0; i < 2; ++i) {
		if (hasLowerBound && (rank % 2 == i)) {
			MPI_Sendrecv_replace(lowerBound, (int)rowSize, MPI_BYTE, rank + 1, 0, rank + 1, 0, MPI_COMM_WORLD,
				&status);
		}
		else if (hasUpperBound & (rank % 2 != i)) {
			MPI_Sendrecv_replace(upperBound, (int)rowSize, MPI_BYTE, rank - 1, 0, rank - 1, 0, MPI_COMM_WORLD,
				&status);
		}
	}

	//��������� ���������������
	for (int row = 0; row < pg->height; row++) {
		for (int column = 0; column < pg->width; column++) {
			fillTempFieldForNode(pg, upperBound, hasUpperBound, lowerBound, hasLowerBound, row, column);
		}
	}
}

/**
* ��������� ���� tmp � ���� �����
* @param grid �����
* @param upperBound ������� �������
* @param hasUpperBound ���� �� ������� ������� � �����
* @param lowerBound ������ �������
* @param hasLowerBound ���� �� ������ ������� � �����
* @param row ������ ����
* @param column ������� ����
*/
void fillTempFieldForNode(const Grid *grid, const GridNode *upperBound, int hasUpperBound, const GridNode *lowerBound,
	int hasLowerBound, int row, int column) {
	GridNode *currentNode = &grid->nodes[row][column];
	double *tmp = currentNode->tmp;
	tmp[0] = currentNode->particleDistribution[0];

	int firstRow = row == 0;
	int firstColumn = column == 0;
	int lastRow = row == grid->height - 1;
	int lastColumn = column == grid->width - 1;

	if (firstRow) {
		if (hasUpperBound) {
			tmp[4] = upperBound[column].particleDistribution[4];
		}
		else {
			tmp[4] = currentNode->particleDistribution[2];
		}
	}
	else {
		tmp[4] = grid->nodes[row - 1][column].particleDistribution[4];
	}
	if (firstColumn) {
		tmp[1] = currentNode->particleDistribution[3];
	}
	else {
		tmp[1] = grid->nodes[row][column - 1].particleDistribution[1];
	}
	if (lastRow) {
		if (hasLowerBound) {
			tmp[2] = lowerBound[column].particleDistribution[2];
		}
		else {
			tmp[2] = currentNode->particleDistribution[4];
		}
	}
	else {
		tmp[2] = grid->nodes[row + 1][column].particleDistribution[2];
	}
	if (lastColumn) {
		tmp[3] = currentNode->particleDistribution[1];
	}
	else {
		tmp[3] = grid->nodes[row][column + 1].particleDistribution[3];
	}
	if (firstRow || firstColumn) {
		if (!firstColumn && hasUpperBound) {
			tmp[8] = upperBound[column - 1].particleDistribution[8];

		}
		else {
			tmp[8] = currentNode->particleDistribution[6];
		}
	}
	else {
		tmp[8] = grid->nodes[row - 1][column - 1].particleDistribution[8];
	}
	if (lastRow || lastColumn) {
		if (!lastColumn && hasLowerBound) {
			tmp[6] = lowerBound[column + 1].particleDistribution[6];
		}
		else {
			tmp[6] = currentNode->particleDistribution[8];
		}
	}
	else {
		tmp[6] = grid->nodes[row + 1][column + 1].particleDistribution[6];
	}
	if (firstRow || lastColumn) {
		if (!lastColumn && hasUpperBound) {
			tmp[7] = upperBound[column + 1].particleDistribution[7];
		}
		else {
			tmp[7] = currentNode->particleDistribution[5];
		}
	}
	else {
		tmp[7] = grid->nodes[row - 1][column + 1].particleDistribution[7];
	}
	if (lastRow || firstColumn) {
		if (!firstColumn && hasLowerBound) {
			tmp[5] = lowerBound[column - 1].particleDistribution[5];
		}
		else {
			tmp[5] = currentNode->particleDistribution[7];
		}
	}
	else {
		tmp[5] = grid->nodes[row + 1][column - 1].particleDistribution[5];
	}
}

/**
* @param tempDistribution �������� ������������� � �����, ���������� �� ����� ���� Streaming
* @param equilibriumDistribution ����������� ������������� �� ������
* @param relaxationTime ����� ���������� ����
* @param result ����� ������������� ������
*/
void updateDistribution(double *tempDistribution,
	double *equilibriumDistribution,
	double relaxationTime,
	double *result) {
	for (int direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
		result[direction] = tempDistribution[direction] +
			(equilibriumDistribution[direction] - tempDistribution[direction]) / relaxationTime;
	}
}

void Collide(Grid *pg) {
	//��������� ������������
	for (int row = 0; row < pg->height; ++row) {
		for (int column = 0; column < pg->width; ++column) {
			GridNode *currentNode = &pg->nodes[row][column];
			// ���������.
			double density = calculateDensity(currentNode->tmp);
			// �������� � �����
			double velocity[2];
			calculateVelocity(currentNode->tmp, density, pg->latticeSpeed, velocity);
			double equilibriumDistribution[LATTICE_DIRECTIONS];
			calculateEquilibriumDistribution(pg->latticeSpeed, density, velocity, equilibriumDistribution);
			// ����� �������������
			updateDistribution(currentNode->tmp, equilibriumDistribution, pg->relaxationTime,
				currentNode->particleDistribution);
		}
	}
}
#pragma endregion
#pragma region Initial data
void initSimulationParameters(const int argc, char **argv,
                              const int worldSize,
                              double *speed, double *relaxationTime, int *totalTime, int *snapshotRate,
                              int *gridWidth) {
    if (argc < 5) {
        printf(
        "usage: boltzman <lattice-speed>  <relaxation-time>  <simulation-time>  <snapshot-rate> (optional <grid-width>)\n");
        exit(1);
    }

    if (sscanf(argv[1], "%lf", speed) != 1) {
        fprintf(stderr, "speed is not double");
        exit(1);
    }
    if (sscanf(argv[2], "%lf", relaxationTime) != 1) {
        fprintf(stderr, "relaxation time is not double");
        exit(1);
    }
    if (sscanf(argv[3], "%i", totalTime) != 1) {
        fprintf(stderr, "simulation time is not integer");
        exit(1);
    }

    if (sscanf(argv[4], "%i", snapshotRate) != 1) {
        fprintf(stderr, "snapshot rate is not integer");
        exit(1);
    }

    if (argc == 6) {
        if (sscanf(argv[5], "%i", gridWidth) != 1) {
            fprintf(stderr, "grid width is not integer");
            exit(1);
        }
    } else {
        // 100�� �� ������ �����������
        *gridWidth = minimumRowCount(sizeof(GridNode), worldSize - 1, 100 * 1024 * 1024);
    }
}

double generateNormalizedRandom() { return rand() / (double)RAND_MAX; }
/**
 * @param from ������, ������� ������������
 * @param to ������, �� ������� ����� �������������
 * @return ���������� ������� ���� ����� ��������� � ����, ��� 0
 */
double tangentProjectionCubed(double *from, double *to) {
    //��� ��� ����� � �������� ������� to ������ ������������ �������,
    //����� ������ �������� ������������ ������ �� ��������
    double cos = cosBetweenVectors(from, to);
    return cos > 0 ? pow(cos, 3) : generateNormalizedRandom() / 20;
}

/**
 * ���������� ������������� ������ �� ������������ � ����� ��� ������������ �������.
 * @param centerOfGrid ����� �������
 * @param row ������
 * @param column �������
 * @param result ������������� ������ �� ������������ � ������ �����.
 */
void generateTwisterData(double *centerOfGrid, int row, int column, double *result) {
    double perpendicular[2];
    perpendicular[0] = centerOfGrid[1] - column;
    perpendicular[1] = row - centerOfGrid[0];
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = tangentProjectionCubed(perpendicular, (double *) elementalVectors[direction]);
    }
}

void InitGrid(Grid *pg, int gridSize, RowBounds bounds, double latticeSpeed, double relaxationTime) {
    double centerOfGrid = (gridSize - 1.) / 2;
    double center[2] = {centerOfGrid, centerOfGrid};
    //������������� �������
    pg->width = gridSize;
    pg->height = bounds.last - bounds.first + 1;
    pg->nodes = calloc((size_t) pg->height, sizeof(GridNode *));
    pg->latticeSpeed = latticeSpeed;
    pg->relaxationTime = relaxationTime;
    int row;
    for (row = 0; row < pg->height; ++row) {
        pg->nodes[row] = calloc((size_t) pg->width, sizeof(GridNode));
        int column;
        for (column = 0; column < pg->width; ++column) {
            GridNode *currentNode = &pg->nodes[row][column];
            generateTwisterData(center, bounds.first + row, column, currentNode->particleDistribution);
        }
    }
}

void FreeGrid(Grid *pg) {
    //������������� �������� �������
}
#pragma endregion

#pragma region snapshots
void getSnapshot(Grid *pg, MacroNode *snapshot) {
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            MacroNode *currentSnapshot = &snapshot[row * pg->width + column];
            GridNode *currentNode = &pg->nodes[row][column];
            currentSnapshot->density = calculateDensity(currentNode->particleDistribution);
            calculateVelocity(currentNode->particleDistribution, currentSnapshot->density, pg->latticeSpeed,
                              currentSnapshot->velocity);
        }
    }
}

void SaveSnapshots(MacroNode *snapshots, int width, int snapshotIndex) {
    //���������� �������
    char fileName[30];
    sprintf(fileName, "snapshot%d.csv", snapshotIndex);
    FILE *file = fopen(fileName, "w");
    fprintf(file, "x,y,Vx,Vy,p\n");
    for (int row = 0; row < width; ++row) {
        for (int column = 0; column < width; ++column) {
            MacroNode *macro = &snapshots[row * width + column];
            fprintf(file, "%d,%d,%f,%f,%f\n", row, column, macro->velocity[0], macro->velocity[1], macro->density);
        }
    }
    fclose(file);
}
#pragma endregion

int main(int argc, char *argv[]) {
    testVectorFunctions();
    testModelFunctions();

    int rank, worldSize, gridWidth, totalTime, snapshotRate;
    double speed, relaxationTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    if (worldSize < 2) {
        printf("world is too small");
        exit(1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    initSimulationParameters(argc, argv, worldSize, &speed, &relaxationTime, &totalTime, &snapshotRate, &gridWidth);

    Grid grid;
    int isMaster = rank == 0;

    //��������� ������ ��� �������� ���������
    int *snapshotSizes = calloc((size_t) worldSize, sizeof(int));
    int *snapshotOffsets = calloc((size_t) worldSize, sizeof(int));
    for (int nonMasterNode = 1; nonMasterNode < worldSize; ++nonMasterNode) {
        RowBounds bounds = getMyBounds(gridWidth, worldSize - 1, nonMasterNode - 1);
        snapshotSizes[nonMasterNode] = (bounds.last - bounds.first + 1) * gridWidth * sizeof(MacroNode);
        snapshotOffsets[nonMasterNode] = bounds.first * gridWidth * sizeof(MacroNode);
    }

    if (!isMaster) {
        RowBounds rowBounds = getMyBounds(gridWidth, worldSize - 1, rank - 1);
        InitGrid(&grid, gridWidth, rowBounds, speed, relaxationTime);
    }
    MacroNode *snapshot;
    if (isMaster) {
        snapshot = calloc((size_t) gridWidth * gridWidth, sizeof(MacroNode));
    } else {
        snapshot = calloc((size_t) grid.height * grid.width, sizeof(MacroNode));
    }
    for (int i = 0; i < totalTime; i++) {
        if (i % snapshotRate == 0) {
            if (!isMaster) {
                getSnapshot(&grid, snapshot);
            }
            MPI_Gatherv(snapshot, isMaster ? 0 : grid.width * grid.height * sizeof(MacroNode), MPI_BYTE, snapshot,
                        snapshotSizes, snapshotOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);
            if (isMaster) {
                SaveSnapshots(snapshot, gridWidth, i / snapshotRate);
            }
        }
        if (!isMaster) {
            Streaming(&grid, rank, worldSize);
            Collide(&grid);
        }
    }
    free(snapshot);
    if (!isMaster) {
        FreeGrid(&grid);
    }
    MPI_Finalize();
}



