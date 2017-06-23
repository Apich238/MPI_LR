#include <mpi.h>
#include <stdlib.h>

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
    double macroscopicDensity;        //макроскопическая плотность
    double macroscopicVelocity[2];        //макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //
    double particleDistribution[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //размеры сетки
    double relaxationTime, latticeSpeed;        //время релаксации и скорость сетки
    GridNode **nodes;
} Grid;

double generateNormalizedRandom() { return rand() / (double) RAND_MAX; }

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

void InitGrid(Grid *pg, int gridSize) {
    //инициализация решётки
    pg->height = pg->width = gridSize;
    pg->nodes = calloc((size_t) gridSize, sizeof(GridNode *));
    pg->latticeSpeed = 1;
    int row;
    for (row = 0; row < pg->height; ++row) {
        pg->nodes[row] = calloc((size_t) gridSize, sizeof(GridNode));
        int column;
        for (column = 0; column < pg->width; ++column) {
            int particleDirection;
            GridNode *currentNode = &pg->nodes[row][column];
            double *probabilitiesOfStreaming = currentNode->particleDistribution;
            double density = 0;
            for (particleDirection = 0; particleDirection < LATTICE_DIRECTIONS; ++particleDirection) {
                probabilitiesOfStreaming[particleDirection] = generateNormalizedRandom();
                density += probabilitiesOfStreaming[particleDirection];
            }
            currentNode->macroscopicDensity = density;
            calculateVelocity(probabilitiesOfStreaming, currentNode->macroscopicDensity, pg->latticeSpeed,
                              currentNode->macroscopicVelocity);
        }
    }
}

void FreeGrid(Grid *pg) {
    //высвобождение ресурсов решётки
}

void Streaming(Grid *pg) {
    //обработка распространения
    int i;
    for (i = 0; i < pg->height; i++) {
        int j;
        for (j = 0; j < pg->width; j++) {

        }
    }
}

void CalcMacro(Grid *pg) {
    //перерасчёт макроскопических плотностей, скоростей, fi eq

}

void Collide(Grid *pg) {
    //обработка столкновений
}

void Snapshot(Grid *pg) {
    //сохранение снимка
}

void SaveSnapshots() {
    //сохранение снимков
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Grid grid;
    int size = 10;
    InitGrid(&grid, 10);
    int totaltime = 10000;
    int snapshoprate = 1000;
    int i;
    for (i = 0; i < totaltime; i++) {
        Streaming(&grid);
        CalcMacro(&grid);
        Collide(&grid);
        if (i % snapshoprate == 0) {
            Snapshot(&grid);
        }
    }
    SaveSnapshots();
    FreeGrid(&grid);
    MPI_Finalize();
}