#include <mpi.h>
#include <stdlib.h>
#include <math.h>

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
    double macroscopicDensity;        //макроскопическа€ плотность
    double macroscopicVelocity[2];        //макроскопическа€ скорость, 0 - горизонтельно, 1 - вертикально
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //
    double particleDistribution[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //размеры сетки
    double relaxationTime, latticeSpeed;        //врем€ релаксации и скорость сетки
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

double tangentProjectionCubed(double *from, double *to) {
    //“ак как здесь в качестве вектора to только элементарные вектора,
    //можно просто умножить элементарный вектор на проекцию
    double cos = cosBetweenVectors(from, to);
    return cos > 0 ? pow(cos,3) : 0;
}

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
    //инициализаци€ решЄтки
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
    //высвобождение ресурсов решЄтки
}

void Streaming(Grid *pg) {
	//обработка распространени€
	//f0 никуда не двигаетс€
	//f1 вправо,f3 влево
	for (int i = 0; i < pg->height; i++) {
		double f = pg->nodes[i][pg->width - 1].particleDistribution[1];
		for (int j = pg->width - 2; j >= 0; j--)
			pg->nodes[i][j + 1].particleDistribution[1] = pg->nodes[i][j].particleDistribution[1];
		pg->nodes[i][0].particleDistribution[1] = pg->nodes[i][0].particleDistribution[3];
		for (int j = 0; j < pg->width - 2; j++)
			pg->nodes[i][j].particleDistribution[3] = pg->nodes[i][j + 1].particleDistribution[3];
		pg->nodes[i][pg->width - 1].particleDistribution[3] = f;
	}
	//f2 вверх, f4 вниз
	for (int j = 0; j < pg->width; j++) {
		double f = pg->nodes[pg->height-1][j].particleDistribution[4];
		for (int i = pg->height - 2; i >= 0; i--)
			pg->nodes[i + 1][j].particleDistribution[4] = pg->nodes[i][j].particleDistribution[4];
		pg->nodes[0][j].particleDistribution[4] = pg->nodes[0][j].particleDistribution[2];
		for (int i = 0; i < pg->height - 2; i++)
			pg->nodes[i][j].particleDistribution[2] = pg->nodes[i + 1][j].particleDistribution[2];
		pg->nodes[pg->height - 1][j].particleDistribution[2] = f;
	}
	//f6 влево вверх, f8 вправо вниз

	//f5 вправо вверх, f7 влево вниз

}

void CalcMacro(Grid *pg) {
    //перерасчЄт макроскопических плотностей, скоростей, fi eq

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
    int speed = 1;
    InitGrid(&grid, size, speed);
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