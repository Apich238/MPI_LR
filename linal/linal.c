#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <stddef.h>
#include <malloc.h>
#include <stdbool.h>
#include <math.h>

//using namespace::std;
int size=1, rank=0;

//структура для элементов списка главных компонент
typedef struct MatrixElem {
	double value;
	int r;
	int c;
} MatrixElem;

typedef struct ColumnHead{
	int id;
	bool enabled;
} ColumnHead;

typedef struct MatrixPart {
	int Rows;
	int Cols;
	bool* RowsEnb;
	ColumnHead* h;
	double** v;
} MatrixPart;

typedef struct IndexInterval {
	int begin;
	int end;
} IndexInterval;

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

//a+=b*c
void AddRow(int len, double* a, double* b, double c) {
	for (int i = 0; i < len; i++)
		a[i] += c*b[i];
}

void WriteResults(double* result, int len) {
	for (int i = 0; i < len; i++) printf("%10.6f ", result[i]);
	printf("\n");
}

void PrintMatr(MatrixPart* m) {
	for (int j = 0; j < m->Cols; j++)
		printf("%10d", m->h[j].id);
	printf("\n");
	for (int i = 0; i < m->Rows; i++) {
		for (int j = 0; j < m->Cols; j++)
			printf("%10.6f", m->v[i][j]);
		printf("\n");
	}
}

void FreeMatrix(MatrixPart* pM) {
	free(pM->h);
	free(pM->RowsEnb);
	for (int i = 0; i < pM->Rows; i++) free(pM->v[i]);
	free(pM->v);
}

double frand(double a,double b) { return a+(b-a)*(rand() / ((double)RAND_MAX)); }

void InitMatrix(MatrixPart * m, int msize, IndexInterval myinterval) {
	m->Rows = msize;
	m->RowsEnb = (bool*)malloc(sizeof(bool)*(m->Rows));
	for (int i = 0; i < m->Rows; i++) m->RowsEnb[i] = true;
	m->v = (double**)malloc(sizeof(double)*(m->Rows));
	m->Cols = myinterval.end - myinterval.begin;
	m->h = (ColumnHead*)malloc((m->Cols) * sizeof(ColumnHead));
	for (int i = 0; i < m->Cols; i++) {
		m->h[i].enabled = true;
		m->h[i].id = myinterval.begin + i;
	}
	if (rank == size - 1) {
		m->h[m->Cols - 1].enabled = false;
		m->h[m->Cols - 1].id = -1;
	}
	for (int i = 0; i < m->Rows; i++)
		m->v[i] = (double*)malloc(sizeof(double)*(m->Cols));
}

void InitValues(MatrixPart* m) {
	srand(rank);
	for (int i = 0; i < m->Rows; i++)
		for (int j = 0; j < m->Cols; j++) {
			m->v[i][j] = frand(-1, 1);
			if (i == m->h[j].id) m->v[i][j] = frand(10, 100);
		}
	if (rank == size - 1)
		for (int i = 0; i < m->Rows; i++)
			m->v[i][m->Cols - 1] = frand(10, 100);
}

void GetMatrix(MatrixPart* res, int sz) {
	IndexInterval MyInterval = IntervalSplit(sz + 1, size, rank);
	InitMatrix(res, sz, MyInterval);
	InitValues(res);
}

void FindLocalMainElement(MatrixPart* pm, MatrixElem* res) {
	res->c = -1;
	res->r = -1;
	res->value = 0;
	for (int i = 0; i < pm->Rows; i++) {
		if (pm->RowsEnb[i])
			for (int j = 0; j < pm->Cols; j++) {
				if (pm->h[j].enabled)
					if (fabs(res->value) < fabs(pm->v[i][j])) {
						res->value = pm->v[i][j];
						res->r = i;
						res->c = pm->h[j].id;
					}
			}
	}
}

void FindMainElement(MatrixPart* m, MatrixElem* res, bool* my, double* mcol) {
	MatrixElem locmax;
	FindLocalMainElement(m, &locmax);
	struct { double val; int rank; } localmax, globalmax;
	localmax.rank = rank;
	if (locmax.c < 0) localmax.val = -1;
	else localmax.val = fabs(m->v[locmax.r][locmax.c - (m->h[0].id)]);
	MPI_Allreduce(&localmax, &globalmax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	*my = (globalmax.rank == rank);
	int index[2] = { locmax.r, locmax.c };
	MPI_Bcast(index, 2, MPI_INT, globalmax.rank, MPI_COMM_WORLD);
	res->r = index[0];
	res->c = index[1];
	double a = locmax.value;
	MPI_Bcast(&a, 1, MPI_DOUBLE, globalmax.rank, MPI_COMM_WORLD);
	res->value = a;
	if (*my) {
		int t = res->c - m->h[0].id;
		for (int i = 0; i < m->Rows; i++)
			mcol[i] = m->v[i][t];
	}
	MPI_Bcast(mcol, m->Rows, MPI_DOUBLE, globalmax.rank, MPI_COMM_WORLD);
}

void GaussForward(MatrixPart* m, MatrixElem* mains) {
	MatrixElem currmax;
	bool MyElem;
	double* mcol = (double*)malloc(m->Rows * sizeof(double));
	for (int i = 0; i < m->Rows; i++) {
		FindMainElement(m, &currmax, &MyElem, mcol);
		mains[i] = currmax;
		if (MyElem) m->h[currmax.c - m->h[0].id].enabled = false;
		m->RowsEnb[currmax.r] = false;
		for (int ri = 0; ri < m->Rows; ri++)
			if (m->RowsEnb[ri])
				AddRow(m->Cols, m->v[ri], m->v[currmax.r], -mcol[ri] / currmax.value);
	}
	free(mcol);
}

void GetCol(int id, double* col, MatrixPart* m) {
	if (id > -1) {
		bool my = m->h[0].id <= id;
		if (rank == size - 1) my = my && id <= m->h[m->Cols - 2].id;
		else my = my&&id <= m->h[m->Cols - 1].id;
		struct { int v; int r; } loc, all;
		if (my) {
			int t = id - m->h[0].id;
			for (int i = 0; i < m->Rows; i++)
				col[i] = m->v[i][t];
			loc.v = 1;
		}
		loc.r = rank;
		MPI_Allreduce(&loc, &all, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
		MPI_Bcast(col, m->Rows, MPI_DOUBLE, all.r, MPI_COMM_WORLD);
	}
	else {
		if (rank == size - 1) {
			for (int i = 0; i < m->Rows; i++)
				col[i] = m->v[i][m->Cols - 1];
		}
		MPI_Bcast(col, m->Rows, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
	}
}

void GaussBack(MatrixPart* m, MatrixElem* mains) {
	double* col = (double*)malloc(sizeof(double)*(m->Rows));
	for (int a = m->Rows - 1; -1 <= a; a--) {
		MatrixElem mel = mains[a];
		GetCol(mel.c, col, m);
		for (int b = 0; b < a; b++)
			AddRow(m->Cols, m->v[mains[b].r], m->v[mel.r], -col[mains[b].r] / mel.value);
	}
	free(col);
}

typedef struct ArrayMember {
	double value;
	int index;
} ArrayMember;

void GetResult(MatrixPart* m, MatrixElem* mains, double* res) {
	double* scol = malloc(sizeof(double)*(m->Rows));
	GetCol(-1, scol, m);
	int cc = m->Cols;
	if (rank == size - 1) cc--;
	ArrayMember* locs = (ArrayMember*)malloc(sizeof(ArrayMember)*cc), *root = NULL;
	for (int i = 0, k = 0; i < m->Rows && k < m->Cols; i++) {
		MatrixElem el = mains[i];
		if (m->h[0].id <= el.c && el.c <= m->h[cc - 1].id) {
			locs[k].index = el.c;
			locs[k].value = scol[el.r] / el.value;
			k++;
		}
	}
	int *szs = NULL, *displs = NULL;
	if (rank == 0) {
		root = (ArrayMember*)malloc(sizeof(ArrayMember)*m->Rows);
		szs = (int*)malloc(sizeof(int)*size);
		displs = (int*)malloc(sizeof(int)*size);
	}
	MPI_Gather(&cc, 1, MPI_INT, szs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		displs[0] = 0;
		for (int i = 1; i < size; i++) displs[i] = displs[i - 1] + szs[i - 1];
	}
	MPI_Gatherv(locs, cc, MPI_DOUBLE_INT, root, szs, displs, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		for (int i = 0; i < m->Rows; i++)
			res[root[i].index] = root[i].value;
		free(root);
		free(displs);
		free(szs);
	}
	free(locs);
	free(scol);
}

void Gauss(MatrixPart* pM, double* solution) {
	MatrixElem* mains = (MatrixElem*)malloc(sizeof(MatrixElem)*(pM->Rows));
	GaussForward(pM, mains);
	GaussBack(pM, mains);
	GetResult(pM, mains, solution);
	free(mains);
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double tCreate = 0, tEvaluate = 0;
	int sz = (int)(10+10240*sqrt(sizeof(double)*1.0/size));
	MatrixPart MyPart;
	double t = MPI_Wtime();
	GetMatrix(&MyPart, sz);
	tCreate = MPI_Wtime() - t;
	//PrintMatr(&MyPart);
	double* sol = NULL;
	if (rank == 0) sol = (double*)malloc(sizeof(double)*(sz));
	t = MPI_Wtime();
	Gauss(&MyPart, sol);
	tEvaluate += MPI_Wtime() - t;
	//if (rank == 0) WriteResults(sol, sz);
	FreeMatrix(&MyPart);
	if (rank == 0) free(sol);
	if (rank == 0) printf("matrix creation time=%f, gauss method time=%f\n", tCreate, tEvaluate);
	MPI_Finalize();
	return 0;
}