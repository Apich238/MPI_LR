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
    double macroscopicDensity;        //макроскопическая плотность
    double macroscopicVelocity[2];        //макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //равновесное распределение
    double particleDistribution[LATTICE_DIRECTIONS];		//распределения частиц по направлениям
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
void equilibriumDistribution(double latticeSpeed, double density, double *velocity, double *result) {
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

void InitGrid(Grid *pg, int gridSize, double latticeSpeed) {
    double centerOfGrid = (gridSize - 1.) / 2;
    double center[2] = {centerOfGrid, centerOfGrid};
    //инициализация решётки
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
    //высвобождение ресурсов решётки
}

void Streaming(Grid *pg) {
	for (int i = 0; i < pg->height; i++)
		for (int j = 0; j < pg->width; j++) {
			pg->nodes[i][j].tmp[0] = pg->nodes[i][j].particleDistribution[0];
			if (i == 0)
				pg->nodes[i][j].tmp[4] = pg->nodes[i][j].particleDistribution[2];
			else
				pg->nodes[i][j].tmp[4] = pg->nodes[i - 1][j].particleDistribution[4];
			if (j == 0)
				pg->nodes[i][j].tmp[1] = pg->nodes[i][j].particleDistribution[3];
			else
				pg->nodes[i][j].tmp[1] = pg->nodes[i][j - 1].particleDistribution[1];
			if (i == pg->height - 1)
				pg->nodes[i][j].tmp[2] = pg->nodes[i][j].particleDistribution[4];
			else
				pg->nodes[i][j].tmp[2] = pg->nodes[i + 1][j].particleDistribution[2];
			if (j == pg->width - 1)
				pg->nodes[i][j].tmp[3] = pg->nodes[i][j].particleDistribution[1];
			else
				pg->nodes[i][j].tmp[3] = pg->nodes[i][j + 1].particleDistribution[3];
			if (i == 0 || j == 0)
				pg->nodes[i][j].tmp[8] = pg->nodes[i][j].particleDistribution[6];
			else
				pg->nodes[i][j].tmp[8] = pg->nodes[i - 1][j - 1].particleDistribution[8];
			if (i == pg->height - 1 || j == pg->width - 1)
				pg->nodes[i][j].tmp[6] = pg->nodes[i][j].particleDistribution[8];
			else
				pg->nodes[i][j].tmp[6] = pg->nodes[i + 1][j + 1].particleDistribution[6];
			if (i == 0 || j == pg->width - 1)
				pg->nodes[i][j].tmp[7] = pg->nodes[i][j].particleDistribution[7];
			else
				pg->nodes[i][j].tmp[7] = pg->nodes[i - 1][j + 1].particleDistribution[7];
			if (i == pg->height - 1 || j == 0)
				pg->nodes[i][j].tmp[5] = pg->nodes[i][j].particleDistribution[7];
			else
				pg->nodes[i][j].tmp[5] = pg->nodes[i + 1][j - 1].particleDistribution[5];
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
    double speed = 0.5;
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
