#include <stdlib.h>
#include <malloc.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>

#define LATTICE_DIRECTIONS 9

#pragma region Typedefs

typedef struct IndexInterval {
	int begin;
	int end;
} IndexInterval;

typedef struct {
	double density;														//макроскопическая плотность
	double velocity[2];													//макроскопическая скорость {Vx,Vy}
} MacroNode;

typedef struct {
	double particleDistribution[LATTICE_DIRECTIONS];        //распределения частиц по направлениям
	double tmp[LATTICE_DIRECTIONS];
} GridNode;

typedef struct {
	int height, width;				//полные размеры сетки
	int locHeight;					//локальная ширина
	IndexInterval LocInds;			//локальные индексы
	GridNode *nodes;				//узлы
	GridNode *TopTmp, *BottomTmp;	//узлы для обмена данными с соседними процессами
} GridPart;							//локальная часть сетки

typedef struct {
	int locheight;
	int width;
	IndexInterval LocInd;
	MacroNode* nodes;
} SnapshotPart;

typedef struct {
	int height;
	int width;
	MacroNode *nodes;
} Snapshot;

#pragma endregion

#pragma region Macro constants and parameters

int rank, worldSize;

bool isMaster;

double weights[] = {4.0 / 9,											//весовые коэффициенты
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

const double e[LATTICE_DIRECTIONS][2] = { {0,  0},						//элементарные вектора
	{1,0},{0,1},{-1,0},{0,-1},
	{1,1},{-1,1},{-1,-1},{1,-1} };

double relaxationTime;		//время релаксации

double latticeSpeed;		//скорость сетки


void InitGlobal() {
	latticeSpeed = 1;
	relaxationTime = 0.7;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	isMaster = (rank == 0);
}

#pragma endregion

#pragma region распределение строк

//example of rank length distribution
//            9
//*  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  7
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
//0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16
//17
//17*7+9=128

int RankLen(int FullLen, int Size, int Rank) { return  FullLen / Size + ((Rank < (FullLen%Size)) ? 1 : 0); }

IndexInterval IntervalSplit(int len, int Size, int Rank) {
	IndexInterval res;
	res.begin = 0;
	for (int i = 0; i < Rank; i++) res.begin += RankLen(len, Size, i);
	res.end = res.begin + RankLen(len, Size, Rank);
	return res;
}

int minimumRowCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
	return (int)ceil((sqrt((minimumSizeOfSystemPerNode * numberOfComputationalNodes)
		/ dataTypeSizeInBytes)));
}

#pragma endregion

#pragma region vecfun
//={0,0}
void VecZero(double *c) {
	for (int i = 0; i < 2; i++) c[i] = 0;
}
//c=a+b
void VecAdd(double *a, double *b, double *c) {
	int i;
	for (i = 0; i < 2; ++i) {
		c[i] = a[i] + b[i];
	}
}
//c=a*b multiply by scale
void VecScale(double *a, double b, double *c) {
	int i;
	for (i = 0; i < 2; ++i) {
		c[i] = a[i] * b;
	}
}
//=|a|
double VecLen(double *a) {
	return sqrt(a[0]*a[0] + a[1]*a[1]);
}
//=(a,b) scalar product
double VecDot(double *a, double *b) {
	double result = 0;
	int i;
	for (i = 0; i < 2; ++i) {
		result += a[i] * b[i];
	}
	return result;
}
//=(a,b)/(|a|*|b|)
double VecCos(double *a, double *b) {
	return VecDot(a, b) / (VecLen(a) * VecLen(b));
}


#pragma endregion

#pragma region функции, связянные с численным методом

/**
 * u = c/ro * summ(fi*ei,i=0..8)
 * @param f распределение частиц по направлениям
 * @param macroscopicDensity микроскопическая плотность в точке
 * @param latticeSpeed скорость сетки
 * @param result микроскопическая скорость в точке
 */
void calculateVelocity(double *f, double macroscopicDensity, double *result) {
	double temp[2];
	VecZero(temp);
	for (int dir = 0; dir < LATTICE_DIRECTIONS; ++dir) {
		VecScale(e[dir], f[dir], temp);
		VecAdd(result, temp, result);
	}
	VecScale(result, latticeSpeed / macroscopicDensity, result);
}

/**
 * ro = summ(fi,i=0..8)
 * @param particleDistribution распределение частиц по направлениям
 * @return микроскопическая плотность в точке
 */
double calculateDensity(double *f) {
    double density = 0;
    for (int dir = 0; dir < LATTICE_DIRECTIONS; ++dir) {
        density += f[dir];
    }
    return density;
}

/**
 * feq[i]=w[i]*(ro+ro*3/c*(e[i],u)+ro*3/2c^2*(3(e[i],u)^2-(u,u))
 * @param latticeSpeed скорость сетки
 * @param density микроскопическая плотность
 * @param velocity микроскопическая скорость
 * @param result равновесное распределение по направлениям (OUT)
 */
void EquilibriumDistribution(double density, double *velocity, double *result) {
	double u2 = VecDot(velocity, velocity);
	double	a = 3 * density / latticeSpeed;
	double	b = a / (2 * latticeSpeed);
	for (int i = 0; i < LATTICE_DIRECTIONS; ++i) {
		double eiu = VecDot(e[i], velocity);
		result[i] = weights[i] * (density + a*eiu + b*(3 * eiu*eiu - u2));
	}
}


/**
* Заполняет поле tmp в ноде сетки
* @param grid сетка
* @param upperBound верхняя граница
* @param hasUpperBound есть ли верхняя граница у сетки
* @param lowerBound нижняя граница
* @param hasLowerBound есть ли нижняя граница у сетки
* @param row строка ноды
* @param column столбец ноды
*/
void fillTempFieldForNode(const GridPart *pg, bool hasUpperBound, bool hasLowerBound, int locrow, int column) {
	GridNode *currentNode = &pg->nodes[locrow*pg->width + column];
	GridNode* upperBound = pg->TopTmp;
	GridNode* lowerBound = pg->BottomTmp;
	double *tmp = currentNode->tmp;
	tmp[0] = currentNode->particleDistribution[0];

	bool firstRow = locrow == 0;
	bool firstColumn = column == 0;
	bool lastColumn = column == pg->width - 1;
	bool lastRow = locrow == pg->locHeight - 1;

	if (firstRow) {
		if (hasUpperBound) {
			tmp[4] = upperBound[column].particleDistribution[4];
		}
		else {
			tmp[4] = currentNode->particleDistribution[2];
		}
	}
	else {
		tmp[4] = pg->nodes[pg->width*(locrow - 1) + column].particleDistribution[4];
	}
	if (firstColumn) {
		tmp[1] = currentNode->particleDistribution[3];
	}
	else {
		tmp[1] = pg->nodes[pg->width*locrow + column - 1].particleDistribution[1];
	}
	if (lastRow) {
		if (hasLowerBound) {
			tmp[2] = lowerBound[column].particleDistribution[2];
		}
		else {
			tmp[2] = currentNode->particleDistribution[4];
		}
	}
	else {
		tmp[2] = pg->nodes[pg->width*(locrow + 1) + column].particleDistribution[2];
	}
	if (lastColumn) {
		tmp[3] = currentNode->particleDistribution[1];
	}
	else {
		tmp[3] = pg->nodes[pg->width*locrow + column + 1].particleDistribution[3];
	}
	if (firstRow || firstColumn) {
		if (!firstColumn && hasUpperBound) {
			tmp[8] = upperBound[column - 1].particleDistribution[8];
		}
		else {
			tmp[8] = currentNode->particleDistribution[6];
		}
	}
	else {
		tmp[8] = pg->nodes[pg->width*(locrow - 1) + column - 1].particleDistribution[8];
	}
	if (lastRow || lastColumn) {
		if (!lastColumn && hasLowerBound) {
			tmp[6] = lowerBound[column + 1].particleDistribution[6];
		}
		else {
			tmp[6] = currentNode->particleDistribution[8];
		}
	}
	else {
		tmp[6] = pg->nodes[pg->width*(locrow + 1) + column + 1].particleDistribution[6];
	}
	if (firstRow || lastColumn) {
		if (!lastColumn && hasUpperBound) {
			tmp[7] = upperBound[column + 1].particleDistribution[7];
		}
		else {
			tmp[7] = currentNode->particleDistribution[5];
		}
	}
	else {
		tmp[7] = pg->nodes[pg->width*(locrow - 1)+column + 1].particleDistribution[7];
	}
	if (lastRow || firstColumn) {
		if (!firstColumn && hasLowerBound) {
			tmp[5] = lowerBound[column - 1].particleDistribution[5];
		}
		else {
			tmp[5] = currentNode->particleDistribution[7];
		}
	}
	else {
		tmp[5] = pg->nodes[pg->width*(locrow + 1) + column - 1].particleDistribution[5];
	}
}

void Streaming(GridPart *pg) {
	int hasUpperBound = rank > 1;
	int hasLowerBound = rank < (worldSize - 1);
	////Обмен смежными строками сетки
	GridNode* temptop, *tempbot;
	if (hasUpperBound) temptop = calloc(pg->width, sizeof(GridNode));
	if (hasLowerBound) tempbot = calloc(pg->width, sizeof(GridNode));
	MPI_Status status;
	if (hasUpperBound) {
		for (int i = 0; i < pg->width; i++) {
			for (int d = 0; d < LATTICE_DIRECTIONS; d++) {
				temptop[i].particleDistribution[d] = pg->nodes[i].particleDistribution[d];
			}
		}
	}
	if (hasLowerBound) {
		for (int i = 0; i < pg->width; i++) {
			for (int d = 0; d < LATTICE_DIRECTIONS; d++) {
				tempbot[i].particleDistribution[d] = pg->nodes[pg->width*(pg->locHeight - 1) + i].particleDistribution[d];
			}
		}
	}
	MPI_Request r1, r2;
	MPI_Status st;
	if (hasUpperBound) {
		MPI_Isend(temptop, sizeof(GridNode)*pg->width, MPI_BYTE, rank - 1, (rank - 1) * 2, MPI_COMM_WORLD, &r1);
	}
	if (hasLowerBound) {
		MPI_Isend(tempbot, sizeof(GridNode)*pg->width, MPI_BYTE, rank + 1, (rank + 1) * 2 + 1, MPI_COMM_WORLD, &r2);
	}
	if (hasUpperBound) {
		MPI_Recv(pg->TopTmp, sizeof(GridNode)*pg->width, MPI_BYTE, rank - 1, rank * 2 + 1, MPI_COMM_WORLD, &st);
	}
	if (hasLowerBound) {
		MPI_Recv(pg->BottomTmp, sizeof(GridNode)*pg->width, MPI_BYTE, rank + 1, rank * 2, MPI_COMM_WORLD, &st);
	}
	if (hasUpperBound) {
		MPI_Wait(&r1, &st);
		free(temptop);
	}
	if (hasLowerBound) {
		MPI_Wait(&r2, &st);
		free(tempbot);
	}
	//обработка распространения
	for (int row = 0; row < pg->locHeight; row++) {
		for (int column = 0; column < pg->width; column++) {
			fillTempFieldForNode(pg, hasUpperBound, hasLowerBound, row, column);
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
	double *equilibriumDistribution, double *result) {
	double at = 1 / relaxationTime;
	for (int direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
		result[direction] = tempDistribution[direction] +
			(equilibriumDistribution[direction] - tempDistribution[direction])*at;
	}
}

void Collide(GridPart *pg) {
	//обработка столкновений
	for (int row = 0; row < pg->locHeight; ++row) {
		for (int column = 0; column < pg->width; ++column) {
			GridNode *currentNode = &(pg->nodes[row*pg->width + column]);
			MacroNode m;
			m.density = calculateDensity(currentNode->tmp);
			// скорость в точке
			calculateVelocity(currentNode->tmp, m.density, m.velocity);
			double equilibriumDistribution[LATTICE_DIRECTIONS];
			EquilibriumDistribution(m.density, m.velocity, equilibriumDistribution);
			// новое распределение
			updateDistribution(currentNode->tmp, equilibriumDistribution, currentNode->particleDistribution);
		}
	}
}
#pragma endregion

#pragma region начальные данные

//random in [a,b]
double frand(double a, double b) {
	return a + (b - a)* rand() / (double)RAND_MAX;
}

/**
* @param from вектор, который проецируется
* @param to векток, на который нужно спроецировать
* @return позитивный косинус угла между векторами в кубе, или 0
*/
double tangentProjectionCubed(double *from, double *to) {
	//Так как здесь в качестве вектора to только элементарные вектора,
	//можно просто умножить элементарный вектор на проекцию
	double cos = VecCos(from, to);
	return cos > 0 ? pow(cos, 3) : frand(0, 1) / 20;
}

/**
 * Генерирует распределение частиц по направлениям в точке для формирования воронки.
 * @param centerOfGrid центр воронки
 * @param row строка
 * @param column столбец
 * @param result распределение частиц по направлениям в данной точке.
 */
void generateTwisterData(double *centerOfGrid, int row, int column, GridNode* node) {
    double perpendicular[2];
    perpendicular[0] = centerOfGrid[1] - column;
    perpendicular[1] = row - centerOfGrid[0];
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
		node->particleDistribution[direction] =  tangentProjectionCubed(perpendicular, e[direction]);
    }
}

#pragma endregion

#pragma region Resource Control
void InitGrid(GridPart *pg, int gridSize,int commsize,int commrank) {
	IndexInterval lind = IntervalSplit(gridSize, commsize, commrank);
    //инициализация структуры
	pg->LocInds = lind;
	pg->width = gridSize;
	pg->height = gridSize;
	pg->locHeight = lind.end - lind.begin ;
	pg->nodes = calloc((size_t)(pg->locHeight * pg->width), sizeof(GridNode));
	pg->TopTmp = calloc((size_t)pg->width, sizeof(GridNode));
	pg->BottomTmp = calloc((size_t)pg->width, sizeof(GridNode));
    //инициализация значений решётки (начальных условий)
    double centerOfGrid = (gridSize - 1.) / 2;
    double center[2] = {centerOfGrid, centerOfGrid};
    for (int row = 0; row < pg->locHeight; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            GridNode *currentNode = &(pg->nodes[row*(pg->width)+column]);
            generateTwisterData(center, lind.begin + row, column, currentNode);
        }
    }
}

void FreeGrid(GridPart *pg) {
    //высвобождение ресурсов решётки
	free(pg->nodes);
	free(pg->TopTmp);
	free(pg->BottomTmp);
}
#pragma endregion

#pragma region Snapshots
int *SnapshotSizes;//размеры в байтах всех частей снэпшотов (у мастера - 0)
int *SnapshotOffsets;//смещения
void InitSnapshotParams(int width,int height) {
	SnapshotSizes = calloc(worldSize, sizeof(int));
	SnapshotOffsets = calloc(worldSize, sizeof(int));
	SnapshotSizes[0] = 0;
	SnapshotOffsets[0] = 0;
	for (int i = 1; i < worldSize; i++) {
		IndexInterval iinds = IntervalSplit(height, worldSize - 1, rank - 1);
		SnapshotSizes[i] = sizeof(MacroNode)* width*(iinds.end - iinds.begin);
		SnapshotOffsets[i] = SnapshotOffsets[i - 1] + SnapshotSizes[i - 1];
	}
}
void FreeParams() {
	free(SnapshotOffsets);
	free(SnapshotSizes);
}
//функции для сборщика
//инициирует 1 снэпшот, в который будут складываться все части
 void InitSnapshot(int height, int width,Snapshot* s) {
	s->height = height;
	s->width = width;
	s->nodes = calloc((size_t)(height*width), sizeof(MacroNode));
}
//освобождает
void FreeSnapshot(Snapshot* ps) {
	free(ps->nodes);
}
//инициирует 1 часть снэпшота
void InitSnapshotPart(int locHeight, int width,SnapshotPart* s) {
	s->locheight = locHeight;
	s->width = width;
	s->nodes = calloc((size_t)(locHeight*width), sizeof(MacroNode));
}
//освобождает
void FreeSnapshotPart(SnapshotPart* ps) {
	free(ps->nodes);
}
//сборка снэпштов
void CollectSnapshot(SnapshotPart *myPart, Snapshot *Res) {
	MacroNode *sbuf=NULL,* rbuf=NULL;
	MacroNode t[1];
	t[0].density = 0;
	t[0].velocity[0] = 0;
	t[0].velocity[1] = 0;
	if (isMaster) {
		sbuf = t;
		rbuf = Res->nodes;
	}
	else {
		sbuf = myPart->nodes;
		rbuf = NULL;
	}
	MPI_Gatherv(sbuf, SnapshotSizes[rank], MPI_BYTE, rbuf, SnapshotSizes, SnapshotOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);
}
void MakeSnapshot(GridPart *pg, SnapshotPart *sp) {
	for (int row = 0; row < pg->locHeight; ++row) {
		for (int column = 0; column < pg->width; ++column) {
			MacroNode *currM = &(sp->nodes[row * pg->width + column]);
			GridNode *currN = &(pg->nodes[row * pg->width + column]);
			currM->density = calculateDensity(currN->particleDistribution);
			calculateVelocity(currN->particleDistribution, currM->density, currM->velocity);
		}
	}
}
void SaveSnapshot(Snapshot* s, int snapshotIndex) {
	//сохранение снимков
	char fileName[30];
	sprintf(fileName, "snapshot%d.csv", snapshotIndex);
	FILE *file = fopen(fileName, "w");
	fprintf(file, "x,y,Vx,Vy,p\n");
	for (int row = 0; row < s->height; ++row) {
		for (int column = 0; column < s->width; ++column) {
			MacroNode *macro = &(s->nodes[row * s->width + column]);
			fprintf(file, "%d,%d,%f,%f,%f\n", row, column, macro->velocity[0], macro->velocity[1], macro->density);
			//printf("%d,%d,%f,%f,%f\n", row, column, macro->velocity[0], macro->velocity[1], macro->density);
		}
	}
	fclose(file);
}
#pragma endregion

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	InitGlobal();
	GridPart grid;
	int gridSize = 3;//= minimumRowCount(sizeof(GridNode), worldSize - 1, 10 * 1024);
	if (!isMaster) InitGrid(&grid, gridSize, worldSize - 1, rank - 1);
	Snapshot s ;
	SnapshotPart sp ;
	InitSnapshotParams(gridSize, gridSize);
	if (isMaster)  InitSnapshot(gridSize, gridSize,&s);
	else  InitSnapshotPart(grid.locHeight, grid.width,&sp);
	int totaltime = 3;
	int snapshotrate = 1;
	for (int i = 0; i < totaltime; i++) {
		if (0 == i%snapshotrate) {
			if (!isMaster) {
				MakeSnapshot(&grid, &sp);
			}
			CollectSnapshot(&sp, &s);
			if (isMaster) {
				SaveSnapshot(&s,i/snapshotrate);
			}
		}
		if (!isMaster) {
			Streaming(&grid);
			Collide(&grid);
		}
	}
	if (isMaster) FreeSnapshot(&s);
	else FreeSnapshotPart(&sp);
	FreeParams();
	if (!isMaster) FreeGrid(&grid);
	MPI_Finalize();
}
