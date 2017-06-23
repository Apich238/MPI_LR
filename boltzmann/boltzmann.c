#include <mpi.h>

int rank, size;

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

typedef struct GridNode {
    double macroscopicDensity;        //���������������� ���������
    double macroscopicVelocity[2];        //���������������� ��������, 0 - �������������, 1 - �����������
    double equilibriumDistribution[9];        //
    double particleDistribution[9];
} GridNode;

typedef struct Grid {
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
    for (int i = 0; i < pg->height; i++) {
        for (int j = 0; j < pg->width; j++) {

        }
    }
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