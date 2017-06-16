#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <stddef.h>
//#include <malloc.h>
#include <math.h>

using namespace::std;
int size=1, rank=0;

//структура для элементов списка главных компонент
struct MatrixIndex {
	int r;
	int c;
};

struct MatrixRow {
	int index;
	bool flag;
	double* values;
};

struct MatrixPart {
	int RowCount, ColumnCount;
	bool* ColumnFlags;
	MatrixRow* Rows;
};

struct IndexInterval {
	int begin;
	int end;
};

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

void PrintMatrix(MatrixPart* pm) {
	for (int i = 0; i < pm->RowCount; i++) {
		printf("R=%4d:", pm->Rows[i].index);
		for (int j = 0; j < pm->ColumnCount-1; j++)
			printf(" %8.4f", pm->Rows[i].values[j]);
		printf(" | %8.4f\n", pm->Rows[i].values[pm->ColumnCount-1]);
	}
}

void freeMatrix(MatrixPart* pM) {
	for (int i = 0; i < pM->RowCount; i++) {
		free(pM->Rows[i].values);
	}
	free(pM->Rows);
	free(pM->ColumnFlags);
}

double frand() { return rand() / ((double)RAND_MAX); }

void GenerateMatrix(MatrixPart * m, int msize, IndexInterval myinterval) {
	m->RowCount = myinterval.end - myinterval.begin;
	m->ColumnCount = msize + 1;
	m->Rows = (MatrixRow*)malloc(sizeof(MatrixRow)*(m->RowCount));
	m->ColumnFlags = (bool*)malloc(sizeof(bool)*(m->ColumnCount - 1));
	srand(rank);
	for (int i = 0, k = myinterval.begin; i < m->RowCount; i++, k++) {
		m->Rows[i].index = k;
		m->Rows[i].values = (double*)malloc(sizeof(double)*(m->ColumnCount));
		for (int j = 0; j < m->ColumnCount; j++) {
			m->Rows[i].values[j] = frand();
			if (m->Rows[i].index == j)
				m->Rows[i].values[j] = 10 * m->Rows[i].values[j] + 30;
			if (j == m->ColumnCount - 1)	m->Rows[i].values[j] = 0;
		}
	}
	srand(0);
	for (int j = 0; j < m->ColumnCount - 1; j++) {
		double c = frand();
		for (int i = 0; i < m->RowCount; i++) m->Rows[i].values[m->ColumnCount - 1] += m->Rows[i].values[j] * c;
	}
}

void ReadInput(MatrixPart* res, int sz) {
	IndexInterval MyInterval= IntervalSplit(sz, size, rank);
	GenerateMatrix(res, sz,MyInterval);
}

void FindLocalMainElement(MatrixPart* pm, MatrixIndex* localmain) {
	MatrixIndex found;
	found.r = -1;
	found.c = -1;
	double lastabs = -1;
	int colmax = pm->ColumnCount - 1;
	for (int i = 0; i < pm->RowCount; i++)
		if (pm->Rows[i].flag)
			for (int j = 0; j < colmax; j++)
				if (pm->ColumnFlags[j])
					if (fabs(pm->Rows[i].values[j]) > lastabs) {
						found.r = i;
						found.c = j;
						lastabs = fabs(pm->Rows[i].values[j]);
					}
	*localmain = found;
	if (found.c >= 0)
		localmain->r += pm->Rows[0].index;
}

void FindMainElement(MatrixPart* pm, MatrixIndex* index, double* values, int *MERank) {
	MatrixIndex locmaxind;
	FindLocalMainElement(pm, &locmaxind);
	struct { double val; int rank; } localmax, globalmax;
	localmax.rank = rank;
	if (locmaxind.c < 0) localmax.val = -1;
	else localmax.val = fabs(pm->Rows[locmaxind.r - pm->Rows[0].index].values[locmaxind.c]);
	MPI_Allreduce(&localmax, &globalmax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	*MERank = globalmax.rank;
	if (rank == globalmax.rank) {
		*index = locmaxind;
		for (int i = 0; i < pm->ColumnCount; i++)
			values[i] = pm->Rows[locmaxind.r - pm->Rows[0].index].values[i];
	}
	MPI_Bcast(index, 2, MPI_INT, globalmax.rank, MPI_COMM_WORLD);
	MPI_Bcast(values, pm->ColumnCount, MPI_DOUBLE, globalmax.rank, MPI_COMM_WORLD);
}

int GetTotalRows(MatrixPart *mp) {
	int loc__totalrows = mp->Rows[mp->RowCount - 1].index + 1, totalrows;
	MPI_Allreduce(&loc__totalrows, &totalrows, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	return totalrows;
}

void GaussForward(MatrixPart* mp, MatrixIndex* mains, int totalrows) {
	double* tmpValues = (double*)malloc(sizeof(double)*mp->ColumnCount);
	MatrixIndex MELInd;
	int MELRank;
	for (int i = 0; i < totalrows; i++) {
		//find main element
		FindMainElement(mp, &MELInd, tmpValues, &MELRank);
		mains[i] = MELInd;
		mp->ColumnFlags[MELInd.c] = false;
		if (rank == MELRank)
			mp->Rows[MELInd.r - mp->Rows[0].index].flag = false;
		//process matrix
		for (int r = 0; r < mp->RowCount; r++)
			if (mp->Rows[r].flag)
				AddRow(mp->ColumnCount, mp->Rows[r].values, tmpValues, -mp->Rows[r].values[MELInd.c] / tmpValues[MELInd.c]);
	}
	free(tmpValues);
}

bool ContainsRow(MatrixPart *m, int r) {
	int a = m->Rows[0].index, b = a + m->RowCount;
	return a <= r && r < b;
}

int RelIndex(MatrixPart* m, int GlobIndex) {
	return GlobIndex - m->Rows[0].index;
}

void GetRow(MatrixPart* m, int index, double* values) {
	struct { int val; int loc; } loc, glob;
	loc.loc = rank;
	loc.val = ContainsRow(m, index) ? 1 : 0;
	MPI_Allreduce(&loc, &glob, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
	int srcrank;
	srcrank = glob.loc;
	if (rank == srcrank)
		for (int i = 0; i < m->ColumnCount; i++)
			values[i] = m->Rows[RelIndex(m, index)].values[i];
	MPI_Bcast(values, m->ColumnCount, MPI_DOUBLE, srcrank, MPI_COMM_WORLD);
}

void GaussBack(MatrixPart* m, MatrixIndex* mains, int TotalRows) {
	double* tmpValues = (double*)malloc(sizeof(double)*m->ColumnCount);
	int a = m->Rows[0].index, b = a + m->RowCount;
	for (int i = TotalRows - 1; i > -1; i--) {
		MatrixIndex mel = mains[i];
		GetRow(m, mel.r, tmpValues);
		for (int j = 0; j < i; j++) {
			MatrixIndex pmel = mains[j];
			if (m->Rows[0].index <= pmel.r && pmel.r < m->Rows[0].index + m->RowCount) {
				int r = pmel.r - m->Rows[0].index;
				AddRow(m->ColumnCount, m->Rows[r].values, tmpValues, -m->Rows[r].values[mel.c] / tmpValues[mel.c]);
			}
		}
	}
	free(tmpValues);
}

struct ArrayMember {
	double value;
	int index;
};

void GetResult(MatrixPart* m, MatrixIndex* mains, int totalrows, double* res) {
	MatrixIndex* locmains = (MatrixIndex*)malloc(sizeof(MatrixIndex) * m->RowCount);
	ArrayMember* locres = (ArrayMember*)malloc(sizeof(ArrayMember) * m->RowCount);
	for (int i = 0, k = 0; i < totalrows; i++)
		if (ContainsRow(m, mains[i].r)) {
			locmains[k] = mains[i];
			k++;
		}
	for (int i = 0; i < m->RowCount; i++) {
		locres[i].index = locmains[i].c;
		int ri = RelIndex(m, locmains[i].r);
		locres[i].value = m->Rows[ri].values[m->ColumnCount - 1] / m->Rows[ri].values[locmains[i].c];
	}
	ArrayMember* allels = 0;
	int *counts = 0, *displs = 0;
	if (rank == 0) {
		allels = (ArrayMember*)malloc(sizeof(ArrayMember)*(m->ColumnCount - 1));
		counts = (int*)malloc(sizeof(int)*size);
		displs = (int*)malloc(sizeof(int)*size);
	}
	MPI_Gather(&(m->RowCount), 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		displs[0] = 0;
		for (int i = 1; i < size; i++)

			displs[i] = displs[i - 1] + counts[i - 1];
	}
	MPI_Gatherv(locres, m->RowCount, MPI_DOUBLE_INT, allels, counts, displs, MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);
	if (rank == 0)
		for (int i = 0; i < m->ColumnCount - 1; i++) {
			res[allels[i].index] = allels[i].value;
		}
	free(locmains);
	free(locres);
	if (rank == 0) {
		free(allels);
		free(counts);
		free(displs);
	}
}

void Gauss(MatrixPart* pMatr, double* solution) {
	int TotalRows = GetTotalRows(pMatr);
	MatrixIndex* mains = (MatrixIndex*)malloc(sizeof(MatrixIndex)*TotalRows);
	GaussForward(pMatr, mains, TotalRows);
	GaussBack(pMatr, mains, TotalRows);
	GetResult(pMatr, mains, TotalRows, solution);
	free(mains);
}

void WriteResults(double* result, int len) {
	for (int i = 0; i < len; i++) printf("%10.6f ", result[i]);
	printf("\n");
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double tCreate = 0, tEvaluate = 0;
	MatrixPart MyPart;
	int reps = 1;
	for (int i = 0; i < reps; i++) {
		double t = MPI_Wtime();
		ReadInput(&MyPart, 1000);
		tCreate += MPI_Wtime() - t;
		double* sol = (double*)malloc(sizeof(double)*(MyPart.ColumnCount - 1));
		t = MPI_Wtime();
		Gauss(&MyPart, sol);
		tEvaluate += MPI_Wtime() - t;
		if (rank == 0) {
			//WriteResults(sol, MyPart.ColumnCount - 1);
		}
		freeMatrix(&MyPart);
		free(sol);
	}
	if (rank == 0) printf("matrix creation time=%f, gauss method time=%f\n", tCreate / reps, tEvaluate / reps);
	MPI_Finalize();
	return 0;
}