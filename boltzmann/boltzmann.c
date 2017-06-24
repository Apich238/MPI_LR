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
    double density;//макроскопическая плотность
    double velocity[2];    //макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
} MacroNode;

typedef struct {
    MacroNode macroParameters;
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //равновесное распределение
    double particleDistribution[LATTICE_DIRECTIONS];        //распределения частиц по направлениям
    double tmp[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //размеры сетки
    double relaxationTime, latticeSpeed;        //время релаксации и скорость сетки
    GridNode **nodes;
} Grid;

void sumVector(double *first, double *second, double *result);

void multiplyVector(double *vector, double multiplier, double *result);

double modulusOfVector(double *vector);

double scalarMultiplication(double *first, double *second);

double cosBetweenVectors(double *first, double *second);

void testVectorFunctions();

void testVectorSum();

void testVectorMultiply();

void testVectorModulus();

void testScalarMultiplication();

/**
 * @param particleDistribution распределение частиц по направлениям
 * @param macroscopicDensity микроскопическая плотность в точке
 * @param latticeSpeed скорость сетки
 * @param result микроскопическая скорость в точке
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
 * @param particleDistribution распределение частиц по направлениям
 * @return микроскопическая плотность в точке
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
 * @param direction направление
 * @param latticeSpeed скорость сетки
 * @param velocity микроскопическая скорость
 * @return Коэффициент для вычисления равновесного распределения по направлениям
 */
double s(int direction, double latticeSpeed, double *velocity) {
    double scalar = scalarMultiplication((double *) elementalVectors[direction], velocity);
    return weights[direction] * 3 *
           (scalar + (3 * pow(scalar, 2) - scalarMultiplication(velocity, velocity)) / (latticeSpeed * 2)) /
           latticeSpeed;
}

/**
 * @param latticeSpeed скорость сетки
 * @param density микроскопическая плотность
 * @param velocity микроскопическая скорость
 * @param result равновесное распределение по направлениям (OUT)
 */
void calculateEquilibriumDistribution(double latticeSpeed, double density, double *velocity, double *result) {
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = (weights[direction] + s(direction, latticeSpeed, velocity)) * density;
    }
}

/**
 * @param from вектор, который проецируется
 * @param to векток, на который нужно спроецировать
 * @return позитивный косинус угла между векторами в кубе, или 0
 */
double tangentProjectionCubed(double *from, double *to) {
    //Так как здесь в качестве вектора to только элементарные вектора,
    //можно просто умножить элементарный вектор на проекцию
    double cos = cosBetweenVectors(from, to);
    return cos > 0 ? pow(cos, 3) : 0;
}

/**
 * Генерирует распределение частиц по направлениям в точке для формирования воронки.
 * @param centerOfGrid центр воронки
 * @param row строка
 * @param column столбец
 * @param result распределение частиц по направлениям в данной точке.
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

void InitGrid(Grid *pg, int gridSize, double latticeSpeed, double relaxationTime) {
    double centerOfGrid = (gridSize - 1.) / 2;
    double center[2] = {centerOfGrid, centerOfGrid};
    //инициализация решётки
    pg->height = pg->width = gridSize;
    pg->nodes = calloc((size_t) gridSize, sizeof(GridNode *));
    pg->latticeSpeed = latticeSpeed;
    pg->relaxationTime = relaxationTime;
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
    //высвобождение ресурсов решётки
}

void Streaming(Grid *pg) {
    //обработка распространения
    for (int row = 0; row < pg->height; row++) {
        for (int column = 0; column < pg->width; column++) {
            GridNode *currentNode = &pg->nodes[row][column];
            double *tmp = currentNode->tmp;
            tmp[0] = currentNode->particleDistribution[0];
            if (row == 0) {
                tmp[4] = currentNode->particleDistribution[2];
            } else {
                tmp[4] = pg->nodes[row - 1][column].particleDistribution[4];
            }
            if (column == 0) {
                tmp[1] = currentNode->particleDistribution[3];
            } else {
                tmp[1] = pg->nodes[row][column - 1].particleDistribution[1];
            }
            if (row == pg->height - 1) {
                tmp[2] = currentNode->particleDistribution[4];
            } else {
                tmp[2] = pg->nodes[row + 1][column].particleDistribution[2];
            }
            if (column == pg->width - 1) {
                tmp[3] = currentNode->particleDistribution[1];
            } else {
                tmp[3] = pg->nodes[row][column + 1].particleDistribution[3];
            }
            if (row == 0 || column == 0) {
                tmp[8] = currentNode->particleDistribution[6];
            } else {
                tmp[8] = pg->nodes[row - 1][column - 1].particleDistribution[8];
            }
            if (row == pg->height - 1 || column == pg->width - 1) {
                tmp[6] = currentNode->particleDistribution[8];
            } else {
                tmp[6] = pg->nodes[row + 1][column + 1].particleDistribution[6];
            }
            if (row == 0 || column == pg->width - 1) {
                tmp[7] = currentNode->particleDistribution[5];
            } else {
                tmp[7] = pg->nodes[row - 1][column + 1].particleDistribution[7];
            }
            if (row == pg->height - 1 || column == 0) {
                tmp[5] = currentNode->particleDistribution[7];
            } else {
                tmp[5] = pg->nodes[row + 1][column - 1].particleDistribution[5];
            }
        }
    }
}

/**
 * @param tempDistribution значение распределения в точке, полученное во время шага Streaming
 * @param equilibriumDistribution равновесное распределение на основе
 * @param relaxationTime время релаксации газа
 * @param result новое распределение частиц
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
    //обработка столкновений
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            GridNode *currentNode = &pg->nodes[row][column];
            // плотность.
            currentNode->macroParameters.density = calculateDensity(currentNode->tmp);
            // скорость в точке
            calculateVelocity(currentNode->tmp, currentNode->macroParameters.density, pg->latticeSpeed,
                              currentNode->macroParameters.velocity);
            double equilibriumDistribution[LATTICE_DIRECTIONS];
            calculateEquilibriumDistribution(pg->latticeSpeed, currentNode->macroParameters.density,
                                             currentNode->macroParameters.velocity, equilibriumDistribution);
            // новое распределение
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
    //сохранение снимков
}

int main(int argc, char *argv[]) {
    testVectorFunctions();
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Grid grid;
    int size = 10;
    double speed = 0.5;
    double relaxationTime = 1.;
    InitGrid(&grid, size, speed, relaxationTime);
    int totalTime = 10000;
    int snapshotRate = 1000;
    int i;
    for (i = 0; i < totalTime; i++) {
        Streaming(&grid);
        Collide(&grid);
        if (i % snapshotRate == 0) {
            MacroNode **snapshot = getSnapshot(&grid);

            //TODO отправить и очистить память для снепшота.
        }
    }
    SaveSnapshots();
    FreeGrid(&grid);
    MPI_Finalize();
}

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

//-----------Тесты-------------------------
void testScalarMultiplication() {
    double first[2] = {1, 2};
    double second[2] = {3, 4};
    double scalar = scalarMultiplication(first, second);
    if (scalar != 1 * 3 + 2 * 4) {
        exit(5);
    }
}

void testVectorModulus() {
    double first[2] = {4, 3};
    double modulus = modulusOfVector(first);
    if (modulus != 5.) {
        exit(4);
    }
}

void testVectorMultiply() {
    double first[2] = {1, 2};
    double result[2];
    multiplyVector(first, 2, result);
    if (result[0] != 2 || result[1] != 4) {
        exit(3);
    }
}

void testVectorSum() {
    double first[2] = {1, 2};
    double second[2] = {3, 4};
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
