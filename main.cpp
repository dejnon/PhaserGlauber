#include <iostream>
#include <algorithm>    // std::swap
#include <math.h>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <cctype>      // old <ctype.h>
#include <sys/types.h>
#include <time.h>
#include <fstream>
#include <unistd.h>
#include <ios>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define SPIN_UP 0
#define SPIN_DOWN 1

#define DISTRIBUTION_WELL 0
#define DISTRIBUTION_GAUSS 1
#define DISTRIBUTION_UNIFORM 2
#define DISTRIBUTION_TRIANGLE 3

// Those below are due to change for different parameters passed
static int LATICE_SIZE = 100;

static double W0 = 0.1;
static int C_MODE = DISTRIBUTION_GAUSS;
static int AVERAGES = 10;
static double C_MIU = 0.5;
static double C_SIGMA = 0.1;

static int MAX_TIME = 1000000;

static bool VERBOSE = false;

clock_t begin, end;
static double time_spent;
static int monte_carlo_steps;


const gsl_rng_type * T;
gsl_rng * r;

double randomGauss(double min=0, double max=1, double sigma=C_SIGMA, double miu=C_MIU) {
    while (true) {
        double new_random = C_MIU+gsl_ran_gaussian_ziggurat(r, C_SIGMA);
        if (min <= new_random && new_random <= max) {
            return new_random;
        }
    }
}

double randomUniform() {
    return gsl_rng_uniform(r);
}


void initAntiferromagnet(short int arr[]) {
    for (int i = 0; i < LATICE_SIZE; i++) {
        arr[i] = (i&1);
    }
}

void initRandom(short int arr[]) {
    for (int i = 0; i < LATICE_SIZE; i++) {
        if (randomUniform()>=0.5) {
            arr[i] = SPIN_UP;
        } else {
            arr[i] = SPIN_DOWN;
        }
    }
}

void initFerromagnet(short int arr[]) {
    for (int i = 0; i < LATICE_SIZE; i++) {
        arr[i] = SPIN_UP;
    }
}

void arrPrint(short int arr[]){
    for (int i = 0; i < LATICE_SIZE; i++) {
        printf(" %d ", arr[i]);
    }
}

short swaped(short val) {
    if (val == SPIN_UP) {
        return SPIN_DOWN;
    } else {
        return SPIN_UP;
    }
}

// Random number from range. Boundaries are inclusive.
int randInRange(int start, int end) {
    return floor(((double)(end-start+1) * randomUniform()) + (double)start);
}

double randInRange(double start, double end) {
    return ((end-start+1.0) * randomUniform()) + start;
}

double distributeWell(double miu, double sigma) {
    return 0.0;
}

double distributeGauss(double miu, double sigma) {
    return randomGauss(sigma=sigma, miu=miu);
}

double distributeUniform(double miu, double sigma) {
    return randInRange(std::max(miu-sigma, 0.0), std::min(miu-sigma, 1.0));
}

double distributeTriangle(double miu, double sigma) {
    double start = std::max(miu-sigma, 0.0);
    double end = std::min(miu+sigma, 1.0);
    double rand = (randomUniform() + randomUniform()) / 2.0;
    return ((end-start) * rand) + start;
}

double bondDensity(short int arr[]) {
    int sum = 0;
    for (int i = 0; i < LATICE_SIZE; i++) {
        int next = (i+1) % LATICE_SIZE;
        sum += 2*abs(arr[i]-arr[next]);
    }
    return (double)sum / (double)(2*LATICE_SIZE);
}

int matnetization(short int arr[]) {
    int sum = 0;
    for (int i = 0; i < LATICE_SIZE; i++) {
        sum += arr[i];
    }
    return sum;
}

double distribution(double miu, double sigma, int mode) {
    switch (mode) {
        case DISTRIBUTION_WELL:
            return distributeWell(miu, sigma);
        case DISTRIBUTION_GAUSS:
            return distributeGauss(miu, sigma);
        case DISTRIBUTION_UNIFORM:
            return distributeUniform(miu, sigma);
        case DISTRIBUTION_TRIANGLE:
            return distributeTriangle(miu, sigma);
        default:
            return 0;
    }
}

std::string distributionName(int mode)
{
    switch (mode) {
        case DISTRIBUTION_WELL:
            return "well";
        case DISTRIBUTION_GAUSS:
            return "gauss";
        case DISTRIBUTION_UNIFORM:
            return "uniform";
        case DISTRIBUTION_TRIANGLE:
            return "triangle";
        default:
            return "-";
    }
}

void updateMonteCarloSteps(bool reset = false) {
    static int spins_updated = 0;
    if (reset) {
        spins_updated = 0;
        monte_carlo_steps = 0;
        return;
    }
    spins_updated++;
    if (spins_updated >= LATICE_SIZE) {
        monte_carlo_steps++;
        spins_updated=0;
    }
}

void displayHelp() {
    std::cout
        << "Usage: \n"
        << "programname [maxt] [l] [w0] [cmean] [csigma] [cmodename]\n"
        << "[maxt]      - Maximal time-step threshold (int)\n"
        << "[l]         - Latice size (int)\n"
        << "[w0]        - Ordering parameter (float), [0,1]\n"
        << "[cmean]     - Mean value of c-parameter, (flt), [0,1]\n"
        << "[csigma]    - C-parameter's standard deviation (or c=[cmean-csigma; cmean+sigma]), (flt), [0,1]\n"
        << "[cmodename] - C-mode: 0-well / 1-gaussian / 2-uniform / 3-triangle, (int), {0,1,2,3}\n"
        << "[verbose]   - Display detailed system progression, (int), {0,1}\n"
        << "[averages]  - How many iterations of given settings, (int)\n"
        << "[file]  - File to write, (str)\n"
    ;
}

int main(int argc, const char * argv[]) {
    std::cout << "==== NEW RUN ====\n";
    if (argc < 9) {
        displayHelp();
        return 0;
    } else if (1 < argc) {
        MAX_TIME    = atoi(argv[1]);
        LATICE_SIZE = atoi(argv[2]);
        W0          = atof(argv[3]);
        C_MIU       = atof(argv[4]);
        C_SIGMA     = atof(argv[5]);
        C_MODE      = atof(argv[6]);
        VERBOSE     = atoi(argv[7]);
        AVERAGES    = atoi(argv[8]);
    } else {
        LATICE_SIZE = 100;
        W0 = 0.1;
        C_MODE = DISTRIBUTION_GAUSS;
        C_MIU = 0.5;
        C_SIGMA = 0.1;
        MAX_TIME = 1000000;
        AVERAGES = 10;
        VERBOSE = false;
    }
    std::ofstream outputFile( argv[9] , std::ios::out | std::ios::app );
    std::cout << "CHECKING FILE\n";
    std::cout << "Writing to: '" << argv[9] << "' \n";
    if (!outputFile) {
        std::cout << argv[9] << "\n";
        printf("ERR");
        outputFile.close();
        return 1;
    }
    
    std::cout << "FILE OK\n";
    
    gsl_rng_env_setup();

    T = gsl_rng_mt19937;
    r = gsl_rng_alloc (T);
    long seed = time (NULL) * getpid();
    gsl_rng_set(r, seed);

    begin = clock();

    std::cout << "MT OK\n";
    
    short * LATICE           = (short *)malloc(LATICE_SIZE*sizeof(short));
    short * NEXT_STEP_LATICE = (short *)malloc(LATICE_SIZE*sizeof(short));

    std::cout << "INIT SIMULATION\n";

    for (int a = 1; a <= AVERAGES; a++) {
        std::cout << "RESET MCS: " << "\n";
        updateMonteCarloSteps(true);
        
        std::cout << "AVERAGE: " << a << "\n";
        initAntiferromagnet(LATICE);
        double start_rho = bondDensity(LATICE);
        double rho=0, sum_rho=0;
        double sum_c = 0;
        int magnetization = 0;
        
        for (int t = 1; t <= MAX_TIME; t++) {
            if (VERBOSE) {
                arrPrint(LATICE);printf("\n");
            }
            rho = bondDensity(LATICE);

            if (rho == 0.0 || t == MAX_TIME) {
                end = clock();
                time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                outputFile << " " << t << " ";
                outputFile << " " << monte_carlo_steps << " ";
                outputFile << " " << matnetization(LATICE) << " ";
                outputFile << " " << rho << "\n";
                std::cout << "cputime:"   << time_spent << "\n";
                break;
            }

            double C = distribution(C_MIU, C_SIGMA, C_MODE);
            sum_c += C;

            int first_i = randInRange(0, LATICE_SIZE);
            int last_i = first_i + (C * LATICE_SIZE);

            for (int i = 0; i < LATICE_SIZE; i++) {

                if (first_i <= i && i <= last_i) {
                    int left  = (i-1) % LATICE_SIZE;
                    int right = (i+1) % LATICE_SIZE;

                    if ( LATICE[left] == LATICE[right] ) {
                        NEXT_STEP_LATICE[i] = LATICE[left];
                        updateMonteCarloSteps();
                    } else if ( W0 > randomUniform() ) {
                        NEXT_STEP_LATICE[i] = swaped(LATICE[i]);
                        updateMonteCarloSteps();
                    }
                } else {
                    NEXT_STEP_LATICE[i] = LATICE[i];
                }

            }

            std::swap(LATICE, NEXT_STEP_LATICE);
        }
    
    }
    gsl_rng_free (r);
    free(LATICE);
    free(NEXT_STEP_LATICE);

    return 0;
}
