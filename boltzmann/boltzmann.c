#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

int rank, size;

double w[] = { 4.0 / 9,
1.0 / 9 ,1.0 / 9 ,1.0 / 9 ,1.0 / 9 ,
1.0 / 36 ,1.0 / 36 ,1.0 / 36 ,1.0 / 36 };

typedef struct GridNode {
	double MDen;
	double MVel[2];
	double feq[9];
	double f[9];
} GridNode;

typedef struct Grid {
	int H, W;
	double tau,c;
	GridNode** nds;
} Grid;

void InitGrid(Grid* pg) {
	//инициализация решётки
}

void FreeGrid(Grid* pg) {
	//высвобождение ресурсов решётки
}

void Streaming(Grid* pg) {
	//обработка распространения
}

void CalcMacro(Grid* pg) {
	//перерасчёт макроскопических плотностей, скоростей, fi eq

}

void Collide(Grid* pg) {
	//обработка столкновений
}

void Snapshot(Grid* pg) {
	//сохранение снимка
}

void SaveSnapshots() {
	//сохранение снимков
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Grid g;
	int size = 10;
	InitGrid(&g, 10);
	int totaltime = 10000;
	int snapshoprate = 1000;
	for (int i = 0; i < totaltime; i++) {
		Streaming(&g);
		CalcMacro(&g);
		Collide(&g);
		if (i%snapshoprate == 0)
			Snapshot(&g);
	}
	SaveSnapshots();
	FreeGrid(&g);
	MPI_Finalize();
}