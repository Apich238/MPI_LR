#include <mpi.h>

int rank, size;

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

typedef struct {
    double macroscopicDensity;        //макроскопическая плотность
    double macroscopicVelocity[2];        //макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
    double equilibriumDistribution[9];        //
    double particleDistribution[9];
} GridNode;

typedef struct {
    int height, width;            //размеры сетки
    double relaxationTime, latticeSpeed;        //время релаксации и скорость сетки
    GridNode **nodes;
} Grid;

void InitGrid(Grid *pg, int gridSize) {
    //инициализация решётки
    pg->height = pg->width = gridSize;
}

void FreeGrid(Grid *pg) {
    //высвобождение ресурсов решётки
}

void Streaming(Grid *pg) {
	//обработка распространения
	//f0 никуда не двигается
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
    for (int i = 0; i < totaltime; i++) {
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