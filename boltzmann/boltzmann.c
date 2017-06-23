#include <mpi.h>

int rank, size;

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

typedef struct {
    double macroscopicDensity;        //���������������� ���������
    double macroscopicVelocity[2];        //���������������� ��������, 0 - �������������, 1 - �����������
    double equilibriumDistribution[9];        //
    double particleDistribution[9];
} GridNode;

typedef struct {
    int height, width;            //������� �����
    double relaxationTime, latticeSpeed;        //����� ���������� � �������� �����
    GridNode **nodes;
} Grid;

void InitGrid(Grid *pg, int gridSize) {
    //������������� �������
    pg->height = pg->width = gridSize;
}

void FreeGrid(Grid *pg) {
    //������������� �������� �������
}

void Streaming(Grid *pg) {
	//��������� ���������������
	//f0 ������ �� ���������
	//f1 ������,f3 �����
	for (int i = 0; i < pg->height; i++) {
		double f = pg->nodes[i][pg->width - 1].particleDistribution[1];
		for (int j = pg->width - 2; j >= 0; j--)
			pg->nodes[i][j + 1].particleDistribution[1] = pg->nodes[i][j].particleDistribution[1];
		pg->nodes[i][0].particleDistribution[1] = pg->nodes[i][0].particleDistribution[3];
		for (int j = 0; j < pg->width - 2; j++)
			pg->nodes[i][j].particleDistribution[3] = pg->nodes[i][j + 1].particleDistribution[3];
		pg->nodes[i][pg->width - 1].particleDistribution[3] = f;
	}
	//f2 �����, f4 ����
	for (int j = 0; j < pg->width; j++) {
		double f = pg->nodes[pg->height-1][j].particleDistribution[4];
		for (int i = pg->height - 2; i >= 0; i--)
			pg->nodes[i + 1][j].particleDistribution[4] = pg->nodes[i][j].particleDistribution[4];
		pg->nodes[0][j].particleDistribution[4] = pg->nodes[0][j].particleDistribution[2];
		for (int i = 0; i < pg->height - 2; i++)
			pg->nodes[i][j].particleDistribution[2] = pg->nodes[i + 1][j].particleDistribution[2];
		pg->nodes[pg->height - 1][j].particleDistribution[2] = f;
	}
	//f6 ����� �����, f8 ������ ����

	//f5 ������ �����, f7 ����� ����

}

void CalcMacro(Grid *pg) {
    //���������� ���������������� ����������, ���������, fi eq

}

void Collide(Grid *pg) {
    //��������� ������������
}

void Snapshot(Grid *pg) {
    //���������� ������
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