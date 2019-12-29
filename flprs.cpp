#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <omp.h>
#include <fstream>

using namespace std;

int numDP = 10000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numCL = 50;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius
int numThreads = 4;

double **demandPoints;  // Geografiniai duomenys
double **distances;
int *X;					// Sprendinys


//=============================================================================

double getTime();
void loadDemandPoints();
void calculateDistances();
void saveToFile();
void loadFromFile();
void randomSolution(int *X);
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);

//=============================================================================

int main() {
    double ts = getTime();          // Algoritmo vykdymo pradzios laikas

    loadDemandPoints();             // Nuskaitomi duomenys
    calculateDistances();           // Skaičiuojami atstumai tarp miestų
    saveToFile();                 // Saugoma i faila
    loadFromFile();                 // Skaitoma is failo

    X = new int[numX];				// Sprndinys
    double u;						// Sprendinio tikslo funkcijos reiksme
    int *bestX = new int[numX];		// Geriausias rastas sprendinys
    double bestU = -1;				// Geriausio sprendinio tikslo funkcijos reiksme

    //----- Pagrindinis ciklas ------------------------------------------------
#pragma omp parallel for num_threads(numThreads)
    for (int iters = 0; iters < 10240; iters++) {
        // Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
        randomSolution(X);
        u = evaluateSolution(X);
        if (u > bestU) {     // Jei geresnis, tai issaugojam kaip geriausia zinoma
            bestU = u;
            for (int i=0; i<numX; i++) bestX[i] = X[i];
        }
    }
    //----- Rezultatu spausdinimas --------------------------------------------

    double tf = getTime();     // Skaiciavimu pabaigos laikas

    cout << "Geriausias sprendinys: ";
    for (int i=0; i<numX; i++) cout << bestX[i] << " ";
    cout << "(" << bestU << ")" << endl << "Skaiciavimo trukme: " << tf-ts << endl;
}

//=============================================================================

void loadDemandPoints() {

    //----- Load demand points ------------------------------------------------
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

void saveToFile() {
    ofstream file;
    file.open("distances.dat");
    for (int i = 0; i < numDP; i++) {
        for (int j = 0; j < numCL; j++) {
            file << distances[i][j];
            if (j != numCL - 1) {
                file << " ";
            }
        }

        if (i != numDP - 1) {
            file << std::endl;
        }
    }
    file.close();
}

void loadFromFile() {
    double ts = getTime();
    distances = new double*[numDP];

    ifstream file;
    file.open("distances.dat");
    for (int i = 0; i < numDP; i++) {
        distances[i] = new double[numCL];
        for (int j = 0; j < numCL; j++) {
            file >> distances[i][j];
        }
    }
    file.close();
    double tf = getTime();
    cout << "Failo skaitymo trukmė: " << tf-ts << endl;
}

void calculateDistances() {
    double ts = getTime();
    distances = new double*[numDP];
#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < numDP; i++) {
        distances[i] = new double[numCL];
#pragma omp parallel for num_threads(numThreads)
        for (int j = 0; j < numCL; j++) {
            distances[i][j] = HaversineDistance(demandPoints[i], demandPoints[j]);
        }
    }
    double tf = getTime();
    cout << "Atstumų skaičiavimo trukmė: " << tf-ts << endl;
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

void randomSolution(int *X) {
    int unique;
    for (int i=0; i<numX; i++) {
        do {
            unique = 1;
            X[i] = (int)((double)rand()/RAND_MAX * numCL);
            for (int j=0; j<i; j++)
                if (X[j] == X[i]) {
                    unique = 0;
                    break;
                }
        } while (unique == 0);
    }
}

//=============================================================================

double evaluateSolution(int *X) {
    double U = 0;
    int bestPF;
    int bestX;
    double d;
    for (int i=0; i<numDP; i++) {
        bestPF = 1e5;
        for (int j=0; j<numPF; j++) {
            d = distances[i][j];
            if (d < bestPF) bestPF = d;
        }
        bestX = 1e5;
        for (int j=0; j<numX; j++) {
            d = distances[i][X[j]];
            if (d < bestX) bestX = d;
        }
        if (bestX < bestPF) U += demandPoints[i][2];
        else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
    }
    return U;
}
#pragma clang diagnostic pop