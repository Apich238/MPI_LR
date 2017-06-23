#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define LATTICE_DIRECTIONS 9

int rank, size;

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
    MacroNode macroParameters;
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //����������� �������������
    double particleDistribution[LATTICE_DIRECTIONS];        //������������� ������ �� ������������
    double tmp[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //������� �����
    double relaxationTime, latticeSpeed;        //����� ���������� � �������� �����
    GridNode **nodes;
} Grid;

void sumVector(double *first, double *second, double *result);

void multiplyVector(double *vector, double multiplier, double *result);

double modulusOfVector(double *vector);

double scalarMultiplication(double *first, double *second);

double cosBetweenVectors(double *first, double *second);

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
        result[direction] = (weights[direction] + s(direction, latticeSpeed, velocity)) * density;
    }
}

/**
 * @param from ������, ������� ������������
 * @param to ������, �� ������� ����� �������������
 * @return ���������� ������� ���� ����� ��������� � ����, ��� 0
 */
double tangentProjectionCubed(double *from, double *to) {
    //��� ��� ����� � �������� ������� to ������ ������������ �������,
    //����� ������ �������� ������������ ������ �� ��������
    double cos = cosBetweenVectors(from, to);
    return cos > 0 ? pow(cos, 3) : 0;
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

void InitGrid(Grid *pg, int gridSize, double latticeSpeed) {
    double centerOfGrid = (gridSize - 1.) / 2;
    double center[2] = {centerOfGrid, centerOfGrid};
    //������������� �������
    pg->height = pg->width = gridSize;
    pg->nodes = calloc((size_t) gridSize, sizeof(GridNode *));
    pg->latticeSpeed = latticeSpeed;
    int row;
    for (row = 0; row < pg->height; ++row) {
        pg->nodes[row] = calloc((size_t) gridSize, sizeof(GridNode));
        int column;
        for (column = 0; column < pg->width; ++column) {
            GridNode *currentNode = &pg->nodes[row][column];
            generateTwisterData(center, row, column, currentNode->particleDistribution);
        }
    }
}

void FreeGrid(Grid *pg) {
    //������������� �������� �������
}

void Streaming(Grid *pg) {
    //��������� ���������������
    for (int i = 0; i < pg->height; i++) {
        for (int j = 0; j < pg->width; j++) {
            pg->nodes[i][j].tmp[0] = pg->nodes[i][j].particleDistribution[0];
            if (i == 0) {
                pg->nodes[i][j].tmp[4] = pg->nodes[i][j].particleDistribution[2];
            } else {
                pg->nodes[i][j].tmp[4] = pg->nodes[i - 1][j].particleDistribution[4];
            }
            if (j == 0) {
                pg->nodes[i][j].tmp[1] = pg->nodes[i][j].particleDistribution[3];
            } else {
                pg->nodes[i][j].tmp[1] = pg->nodes[i][j - 1].particleDistribution[1];
            }
            if (i == pg->height - 1) {
                pg->nodes[i][j].tmp[2] = pg->nodes[i][j].particleDistribution[4];
            } else {
                pg->nodes[i][j].tmp[2] = pg->nodes[i + 1][j].particleDistribution[2];
            }
            if (j == pg->width - 1) {
                pg->nodes[i][j].tmp[3] = pg->nodes[i][j].particleDistribution[1];
            } else {
                pg->nodes[i][j].tmp[3] = pg->nodes[i][j + 1].particleDistribution[3];
            }
            if (i == 0 || j == 0) {
                pg->nodes[i][j].tmp[8] = pg->nodes[i][j].particleDistribution[6];
            } else {
                pg->nodes[i][j].tmp[8] = pg->nodes[i - 1][j - 1].particleDistribution[8];
            }
            if (i == pg->height - 1 || j == pg->width - 1) {
                pg->nodes[i][j].tmp[6] = pg->nodes[i][j].particleDistribution[8];
            } else {
                pg->nodes[i][j].tmp[6] = pg->nodes[i + 1][j + 1].particleDistribution[6];
            }
            if (i == 0 || j == pg->width - 1) {
                pg->nodes[i][j].tmp[7] = pg->nodes[i][j].particleDistribution[5];
            } else {
                pg->nodes[i][j].tmp[7] = pg->nodes[i - 1][j + 1].particleDistribution[7];
            }
            if (i == pg->height - 1 || j == 0) {
                pg->nodes[i][j].tmp[5] = pg->nodes[i][j].particleDistribution[7];
            } else {
                pg->nodes[i][j].tmp[5] = pg->nodes[i + 1][j - 1].particleDistribution[5];
            }
        }
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
            currentNode->macroParameters.density = calculateDensity(currentNode->tmp);
            // �������� � �����
            calculateVelocity(currentNode->tmp, currentNode->macroParameters.density, pg->latticeSpeed,
                              currentNode->macroParameters.velocity);
            double equilibriumDistribution[LATTICE_DIRECTIONS];
            calculateEquilibriumDistribution(pg->latticeSpeed, currentNode->macroParameters.density,
                                             currentNode->macroParameters.velocity, equilibriumDistribution);
            // ����� �������������
            updateDistribution(currentNode->tmp, equilibriumDistribution, pg->relaxationTime,
                               currentNode->particleDistribution);
        }
    }
}

MacroNode **getSnapshot(Grid *pg) {
    MacroNode **snapshot = calloc((size_t) pg->height, sizeof(MacroNode));
    for (int row = 0; row < pg->height; ++row) {
        snapshot[row] = calloc((size_t) pg->width, sizeof(MacroNode));
        for (int column = 0; column < pg->width; ++column) {
            MacroNode *currentSnapshot = &snapshot[row][column];
            GridNode *currentNode = &pg->nodes[row][column];
            memcpy(currentSnapshot, &currentNode->macroParameters, sizeof(MacroNode));
        }
    }
    return snapshot;
}

void SaveSnapshots() {
    //���������� �������
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Grid grid;
    int size = 10;
    double speed = 0.5;
    InitGrid(&grid, size, speed);
    int totalTime = 10000;
    int snapshotRate = 1000;
    int i;
    for (i = 0; i < totalTime; i++) {
        Streaming(&grid);
        Collide(&grid);
        if (i % snapshotRate == 0) {
            MacroNode **snapshot = getSnapshot(&grid);

            //TODO ��������� � �������� ������ ��� ��������.
        }
    }
    SaveSnapshots();
    FreeGrid(&grid);
    MPI_Finalize();
}

void sumVector(double *first, double *second, double *result) {
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = first[0] + second[0];
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
