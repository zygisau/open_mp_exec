#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <sstream>

using namespace std;

int demandPointsCount = 10000;      // Vietoviu skaicius (demand points, max 10000)
int preexistingPointsCount = 5;     // Esanciu objektu skaicius (preexisting facilities)
int candidateLocationsCount = 50;   // Kandidatu naujiems objektams skaicius (candidate locations)
int candidatesCount = 3;            // Nauju objektu skaicius

double **demandPoints;              // Geografiniai duomenys
double **citiesMatrix = nullptr;    // Miestų matrica
int *currentPoint;                  // Sprendinys

//=============================================================================

double getTime();
void loadDemandPoints();
void randomSolution(int *X);
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);
void calculateDistanceMatrix();
void assignDistance(int i, int j);
void initializeMatrix();
void writeToFile(string& fileName);
void writeDistancesToFile(string& fileName);
void ensureCitiesMatrixInitialized();
void readDistancesFromFile(string& fileName);
void readFromFile(string& fileName);

//=============================================================================

int main() {

double ts = getTime();                          // Algoritmo vykdymo pradzios laikas

    loadDemandPoints();                         // Nuskaitomi duomenys

    currentPoint = new int[candidatesCount];    // Sprendinys
    double u;						            // Sprendinio tikslo funkcijos reiksme
    int *bestX = new int[candidatesCount];		// Geriausias rastas sprendinys
    double bestU = -1;	            			// Geriausio sprendinio tikslo funkcijos reiksme

    int nIterations = 10048;                    // Iteraciju skaicius
    int NUM_THREADS = 4;                        // Giju skaicius

    string fileName = "distances_between_cities.txt";
    // writeDistancesToFile(fileName);
    readDistancesFromFile(fileName);
    //----- Pagrindinis ciklas ------------------------------------------------
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        initializeMatrix();
        calculateDistanceMatrix();
        #pragma omp for schedule(static)
        for (int iters=0; iters<nIterations; iters++) {
            // Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
            randomSolution(currentPoint);
            u = evaluateSolution(currentPoint);
            if (u > bestU) {     // Jei geresnis, tai issaugojam kaip geriausia zinoma
                bestU = u;
                for (int i=0; i < candidatesCount; i++) bestX[i] = currentPoint[i];
            }
        }
    }

    //----- Rezultatu spausdinimas --------------------------------------------

    double tf = getTime();     // Skaiciavimu pabaigos laikas

    cout << "Geriausias sprendinys: ";
    for (int i=0; i < candidatesCount; i++) cout << bestX[i] << " ";
    cout << "(" << bestU << ")" << endl << "Skaiciavimo trukme: " << tf-ts << endl;
}

//=============================================================================

void readFromFile(string& fileName) {
    ifstream file;
    file.open(fileName, ofstream::in);
    string line;
    int row = 0, column = 0, value = 0;
    while(std::getline(file, line)) {
        std::stringstream lineStream(line);

        while (lineStream >> value) {
            *(*(citiesMatrix + row) + column) = value;
            column++;
        }
        row++;
        column = 0;
    }
    file.close();
}

//=============================================================================

void readDistancesFromFile(string& fileName) {
    if (citiesMatrix == nullptr) {
        initializeMatrix();
    }
    readFromFile(fileName);
}

//=============================================================================

void writeToFile(string& fileName) {
    ofstream file;
    file.open(fileName, ofstream::out);
    for (int i = 0; i < demandPointsCount; ++i) {
        for (int j = 0; j < demandPointsCount; ++j) {
            file << citiesMatrix[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

//=============================================================================

void ensureCitiesMatrixInitialized() {
    if (citiesMatrix == nullptr) {
        initializeMatrix();
        calculateDistanceMatrix();
    }
}

//=============================================================================

void writeDistancesToFile(string& fileName) {
    ensureCitiesMatrixInitialized();
    writeToFile(fileName);
}

//=============================================================================

void initializeMatrix() {
    citiesMatrix = new double*[demandPointsCount];
    for(int i = 0; i < demandPointsCount; ++i)
        citiesMatrix[i] = new double[demandPointsCount]();
}

//=============================================================================

void assignDistance(int i, int j) {
    if (i == j) {
        citiesMatrix[i][j] = 0;
    } else {
        citiesMatrix[i][j] = HaversineDistance(demandPoints[i], demandPoints[j]);
    }
}

//=============================================================================

void calculateDistanceMatrix() {
    #pragma omp for schedule(static)
    for (int i = 0; i < demandPointsCount; i++) {
        for (int j = 0; j < demandPointsCount; j++) {
            assignDistance(i, j);
        }
    }
}

//=============================================================================

void loadDemandPoints() {

    //----- Load demand points ------------------------------------------------
    FILE *f;
    f = fopen("demandPoints.dat", "r");
    demandPoints = new double*[demandPointsCount];
    for (int i=0; i < demandPointsCount; i++) {
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

void randomSolution(int *X) {
    int unique;
    for (int i=0; i < candidatesCount; i++) {
        do {
            unique = 1;
            X[i] = (int)((double)rand() / RAND_MAX * candidateLocationsCount);
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
    for (int i=0; i < demandPointsCount; i++) {
        bestPF = 1e5;
        for (int j=0; j < preexistingPointsCount; j++) {
            d = citiesMatrix[i][j];
            if (d < bestPF) bestPF = d;
        }
        bestX = 1e5;
        for (int j=0; j < candidatesCount; j++) {
            d = citiesMatrix[i][X[j]];
            if (d < bestX) bestX = d;
        }
        if (bestX < bestPF) U += demandPoints[i][2];
        else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
    }
    return U;
}