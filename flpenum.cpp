#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#include <mpi.h>

using namespace std;

int numDP = 1000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;        // Esanciu objektu skaicius (preexisting facilities)
int numCL = 26;        // Kandidatu naujiems objektams skaicius (candidate locations)

double **demandPoints; // Geografiniai duomenys


//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X, int numX);
int increaseX(int *X, int index, int maxindex, int numX);
int calculateCombinationsCount(int numX);
void calculateValues(int numX, int rank, int comm_size, double ts);

//=============================================================================

int main(int argc, char *argv[]) {
	double ts = getTime();
	ts = getTime();
	int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	loadDemandPoints();

	for (int numX = 3; numX <= 5; numX++) {
		calculateValues(numX, rank, comm_size, ts);
		MPI_Barrier(MPI_COMM_WORLD);
	}

    MPI_Finalize();
	return 0;
}

void calculateValues(int numX, int rank, int comm_size, double ts) {
	int combinationCount = calculateCombinationsCount(numX);
	int *X = new int[numX];
	int bestX[numX];
	for (int i=0; i<numX; i++) {
		X[i] = i;
		bestX[i] = i;
	}
	double u = evaluateSolution(X, numX);
	double bestU = u;
	int r;

	//----- Pagrindinis ciklas ------------------------------------------------
	int i = 0;

	while (true) {
		i++;
		if (increaseX(X, numX-1, numCL, numX)) {
			if (i < (combinationCount / comm_size) * (rank + 1) && i > (combinationCount / comm_size) * (rank)) {
					u = evaluateSolution(X, numX);
					if (u > bestU) {
						bestU = u;
						for (int i=0; i<numX; i++) bestX[i] = X[i];
					}
			} else if(i >= (combinationCount / comm_size) * (rank + 1)) {
				break;
			}
		}
		else break;
	}
	//----- Rezultatu spausdinimas --------------------------------------------

	MPI_Barrier(MPI_COMM_WORLD);

    double *receivedBestU = (double *)malloc(sizeof(double) * comm_size);
	int receivedBestX[comm_size][numX];

    MPI_Gather(&bestU, 1, MPI_DOUBLE, receivedBestU, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&bestX, numX, MPI_INT, &receivedBestX, numX, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		int bU = 0;
		int bI = 0;

		for (int i = 0; i < comm_size; i++) {
			if (receivedBestU[i] > bU) {
				bI = i;
				bU = receivedBestU[i];
			}
		}
	
		double tf = getTime();     // Skaiciavimu pabaigos laikas
		cout << "numX: " << numX << endl;
		cout << "Geriausias sprendinys: ";
		for (int i=0; i<numX; i++) cout << receivedBestX[bI][i] << " ";
		cout << "(" << bU << ")" << endl;
		cout << "Skaiciavimo trukme: " << tf-ts << endl;
	}
}

//-

int calculateCombinationsCount(int numX) {
	int count = numCL;
	int count2 = numX;
	for (int i = 1; i < numX; i++) {
		count *= numCL - i;
		count2 *= numX - i;
	}

	return count / count2;
}

//=============================================================================

void loadDemandPoints() {
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double* a, double* b) {
   double dlon = fabs(a[0] - b[0]);
   double dlat = fabs(a[1] - b[1]);
   double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
   double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
   double d = 6371 * c; 
   return d;
}

//=============================================================================

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

//=============================================================================

double evaluateSolution(int *X, int numX) {
	double U = 0;
	int bestPF;
	int bestX;
	double d;
	for (int i=0; i<numDP; i++) {
		bestPF = 1e5;
		for (int j=0; j<numPF; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF) bestPF = d;
		}
		bestX = 1e5;
		for (int j=0; j<numX; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX) bestX = d;
		}
		if (bestX < bestPF) U += demandPoints[i][2];
		else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
	}
	return U;
}

//=============================================================================

int increaseX(int *X, int index, int maxindex, int numX) {
	if (X[index]+1 < maxindex-(numX-index-1)) {
		X[index]++;
	}
	else {		 
		if ((index == 0) && (X[index]+1 == maxindex-(numX-index-1))) {
			return 0;
		}
		else {
			if (increaseX(X, index-1, maxindex, numX)) X[index] = X[index-1]+1;
			else return 0;
		}	
	}
	return 1;
}