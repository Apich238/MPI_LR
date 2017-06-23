#include <mpi.h>
#include <stdlib.h>

#define LATTICE_DIRECTIONS 9

int rank, size;

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

typedef struct {
    double macroscopicDensity;        //���������������� ���������
    double macroscopicVelocity[2];        //���������������� ��������, 0 - �������������, 1 - �����������
    double equilibriumDistribution[LATTICE_DIRECTIONS];        //
    double particleDistribution[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
    int height, width;            //������� �����
    double relaxationTime, latticeSpeed;        //����� ���������� � �������� �����
    GridNode **nodes;
} Grid;

double generateNormalizedRandom() { return rand() / (double) RAND_MAX; }

void InitGrid(Grid *pg, int gridSize) {
    //������������� �������
    pg->height = pg->width = gridSize;
    pg->nodes = calloc((size_t) gridSize, sizeof(GridNode *));
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
                probabilitiesOfStreaming[particleDirection + 1] = generateNormalizedRandom();
                density += probabilitiesOfStreaming[particleDirection];
            }
            currentNode->macroscopicDensity = density;

        }
    }
}

void FreeGrid(Grid *pg) {
    //������������� �������� �������
}

void Streaming(Grid *pg) {
    //��������� ���������������
    int i;
    for (i = 0; i < pg->height; i++) {
        int j;
        for (j = 0; j < pg->width; j++) {

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