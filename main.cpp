#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

using std::cout, std::endl, std::string, std::ifstream, std::ofstream, std::vector;

int demandPointsCount = 10000;      // Vietoviu skaicius (demand points, max 10000)
int preexistingPointsCount = 5;     // Esanciu objektu skaicius (preexisting facilities)
int candidateLocationsCount = 50;   // Kandidatu naujiems objektams skaicius (candidate locations)
int candidatesCount = 3;            // Nauju objektu skaicius

vector<vector<double>> demandPointsVector;
vector<vector<double>>::iterator demandPoints;              // Geografiniai duomenys
double **citiesMatrix = nullptr;    // Miestų matrica

//=============================================================================

double getTime();
void loadDemandPoints();
void randomSolution(vector<int>& X);
double HaversineDistance(vector<vector<double>>::iterator a, vector<vector<double>>::iterator b);
double evaluateSolution(vector<int>& X);
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

    double ts = getTime();                      // Algoritmo vykdymo pradzios laikas

    loadDemandPoints();                         // Nuskaitomi duomenys

    vector<int> bestX = vector<int>(candidatesCount);		// Geriausias rastas sprendinys
    double bestU = -1;	            			// Geriausio sprendinio tikslo funkcijos reiksme

    int nIterations = 10048;                    // Iteraciju skaicius
    int NUM_THREADS = 4;                        // Giju skaicius

    // double matrixLoadFromFileStart = getTime();
    // string fileName = "distances_between_cities.txt";
    // writeDistancesToFile(fileName);
    // readDistancesFromFile(fileName);
    // double matrixLoadFromFileEnd = getTime();

    omp_set_num_threads(NUM_THREADS);
    
    double matrixCalcStart = getTime();
    initializeMatrix();
    calculateDistanceMatrix();
    double matrixCalcEnd = getTime();
    
    //----- Pagrindinis ciklas ------------------------------------------------
    double mainCalcStart = getTime();
    #pragma omp parallel shared(candidatesCount, nIterations, bestX, bestU) default(none)
    {
        vector<int> threadCurrentPoint(candidatesCount); // Sprendinys
        double threadU;                                     // Sprendinio tikslo funkcijos reiksme
        double threadBestU = -1;

        #pragma omp for schedule(static)
        for (int iters=0; iters<nIterations; iters++) {
            // Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
            randomSolution(threadCurrentPoint);
            threadU = evaluateSolution(threadCurrentPoint);
            if (threadU > threadBestU) {     // Jei geresnis, tai issaugojam kaip geriausia zinoma
                threadBestU = threadU;
            }
        }

        #pragma omp critical
        {
            if (threadBestU > bestU) {
                bestU = threadBestU;
                std::copy(threadCurrentPoint.begin(), threadCurrentPoint.end(), bestX.begin());
            }
        }

    }
    double mainCalcStop = getTime();

    //----- Rezultatu spausdinimas --------------------------------------------

    double tf = getTime();     // Skaiciavimu pabaigos laikas

    cout << "Geriausias sprendinys: ";
    for (int i=0; i < candidatesCount; i++) cout << bestX[i] << " ";
    cout << "(" << bestU << ")" << endl;
    cout << "Matricos skaiciavimo trukme: " << matrixCalcEnd-matrixCalcStart << endl;
    // cout << "Matricos ikelimo is failo trukme: " << matrixLoadFromFileEnd-matrixLoadFromFileStart << endl;
    cout << "Skaiciavimo trukme: " << mainCalcStop-mainCalcStart << endl;
    cout << "Programos veikimo trukme: " << tf-ts << endl;
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
        for (int j = 0; j < candidateLocationsCount; ++j) {
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
        citiesMatrix[i] = new double[candidateLocationsCount]();
}

//=============================================================================

void assignDistance(int i, int j) {
    if (i == j) {
        citiesMatrix[i][j] = 0;
    } else {
        citiesMatrix[i][j] = HaversineDistance(demandPoints+i, demandPoints+j);
    }
}

//=============================================================================

void calculateDistanceMatrix() {
    #pragma omp parallel shared(demandPointsCount, candidateLocationsCount) default(none)
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < demandPointsCount; i++) {
            for (int j = 0; j < candidateLocationsCount; j++) {
                assignDistance(i, j);
            }
        }
    }
}

//=============================================================================

void loadDemandPoints() {

    //----- Load demand points ------------------------------------------------
    ifstream file;
    file.open("demandPoints.dat", ofstream::in);
    string line, substring;
    int column = 0;
    demandPointsVector.reserve(demandPointsCount);
    demandPoints = demandPointsVector.begin();
    for (int i=0; i < demandPointsCount; i++) {
        std::getline(file, line);
        std::stringstream lineStream(line);

        demandPointsVector.emplace_back(3, 0);
        while (std::getline(lineStream, substring, '\t')) {
            demandPointsVector[i][column] = std::stod(substring);
            column++;
        }
        column = 0;
    }
    file.close();
}

//=============================================================================

double HaversineDistance(vector<vector<double>>::iterator a, vector<vector<double>>::iterator b) {
    vector<double> pointA = (*a);
    vector<double> pointB = (*b);
    double dlon = fabs(pointA[0] - pointB[0]);
    double dlat = fabs(pointA[1] - pointB[1]);
    double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(pointA[0]*0.01745) * cos(pointB[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
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

void randomSolution(vector<int>& X) {
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<> random(0, candidateLocationsCount);
    bool unique;
    for (int i=0; i < candidatesCount; i++) {
        do {
            X[i] = random(rng);
            auto it = std::find(X.begin(), (X.begin() + i), X[i]);
            unique = it == (X.begin() + i);
        } while (!unique);
    }
}

//=============================================================================

double evaluateSolution(vector<int>& X) {
    double U = 0;
    int bestPF;
    int bestX;
    double d;
    for (int i=0; i < demandPointsCount; i++) {
        bestPF = 1e5;
        for (int j=0; j < preexistingPointsCount; j++) {
            d = citiesMatrix[i][j];   
            // d = HaversineDistance(demandPoints[i], demandPoints[j]);
            if (d < bestPF) bestPF = d;
        }
        bestX = 1e5;
        for (int j=0; j < candidatesCount; j++) {
            d = citiesMatrix[i][X[j]];
            // d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
            if (d < bestX) bestX = d;
        }
        if (bestX < bestPF) U += demandPointsVector[i][2];
        else if (bestX == bestPF) U += 0.3 * demandPointsVector[i][2];
    }
    return U;
}