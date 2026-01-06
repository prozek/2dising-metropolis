#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <csignal>

/*
* Choose a random double from range [0,1]
*/
double drand01() {
    return double( rand() )/ double( RAND_MAX );
}

class simulation {
public:
    simulation();
    simulation(int inL, double coupl, double field);
    ~simulation();

    double mag();
    double mag2();
    double ener();
    double ener2();
    
    int* genNeigh(int i);
    void getNeighbors(int i, int* Neigh);
    
    void init();
    void MCstep(int site, double b);
    void run(int steps, double b);
    void tempsweep(int steps, double bs, double bf, int bres);
    void writeLat();
    void writeToShm(char* shm_ptr, long step_count);
    int getL() { return L; }
    int getN() { return N; }
     
private:
    int L;
    int N;
    int NN;
    int* lat;
    double J;
    double h;
};

simulation::simulation() {}

simulation::simulation(int inL, double coupl, double field) {
    L = inL;
    N = L*L;
    NN = 4;
    lat = new int[N];
    J = coupl;
    h = field;
}

simulation::~simulation() {}

/*
* Generate positions of nearest neighbors on given site, using time-expensive method.
* Site i is at row = i/L, col = i%L in the LxL lattice
*/
int* simulation::genNeigh(int i) {
    int* Neigh = new int[NN];
    int row = i / L;
    int col = i % L;
    
    // Right neighbor (col+1, with periodic boundary)
    int right_col = (col + 1) % L;
    Neigh[0] = row * L + right_col;
    
    // Left neighbor (col-1, with periodic boundary)
    int left_col = (col - 1 + L) % L;
    Neigh[1] = row * L + left_col;
    
    // Down neighbor (row+1, with periodic boundary)
    int down_row = (row + 1) % L;
    Neigh[2] = down_row * L + col;
    
    // Up neighbor (row-1, with periodic boundary)
    int up_row = (row - 1 + L) % L;
    Neigh[3] = up_row * L + col;
    
    return Neigh;
}

/*
* Optimized inline neighbor lookup - no allocation
* Site i is at row = i/L, col = i%L in the LxL lattice
*/
void simulation::getNeighbors(int i, int* Neigh) {
    int row = i / L;
    int col = i % L;
    
    // Right neighbor (col+1, with periodic boundary)
    int right_col = (col + 1) % L;
    Neigh[0] = row * L + right_col;
    
    // Left neighbor (col-1, with periodic boundary)
    int left_col = (col - 1 + L) % L;
    Neigh[1] = row * L + left_col;
    
    // Down neighbor (row+1, with periodic boundary)
    int down_row = (row + 1) % L;
    Neigh[2] = down_row * L + col;
    
    // Up neighbor (row-1, with periodic boundary)
    int up_row = (row - 1 + L) % L;
    Neigh[3] = up_row * L + col;
}


double simulation::mag() {
	double res = 0;
	for(int i=0;i<N;i++) {
		res += (double) lat[i];	
    }
    return res;
}


double simulation::mag2() {
	double res = 0;
	for(int i=0;i<N;i++) {
	for(int j=0;j<N;j++) {
		res += (double) lat[i]*lat[j]; 
	}}	
	return res;
}


double simulation::ener() {
	double res = 0;
	int Neigh[4];
	for(int i=0;i<N;i++) {			
		getNeighbors(i, Neigh);
		for(int j=0;j<NN;j++) {
			res += (double) lat[i]*lat[Neigh[j]]; 
		}
	}
	return res; 
}


double simulation::ener2() {
	double res = 0;	
	for(int i1=0;i1<N;i1++) {
	for(int j1=0;j1<N;j1++) {
		for(int k1=0;k1<NN;k1++) { 
		for(int l1=0;l1<NN;l1++) {	
			res += (double) lat[i1]*lat[genNeigh(i1)[k1]] * lat[j1]*lat[genNeigh(j1)[l1]];
			}}
		}}
	return res;
}

/*
* Initialize the system in random configuration
*/
void simulation::init() {
    for (int i=0;i<N;i++) {
        double rd = drand01();
        if (rd < 0.5) { lat[i] = 1; } else { lat[i] = -1; }
    }
}

void simulation::MCstep(int i, double b) {
    double deltaE = 0;
    double w;
    int Neigh[4];
    getNeighbors(i, Neigh);
    
    // Coupling term: interaction with 4 nearest neighbors
    for (int j=0;j<NN;j++) { deltaE += 2. * J * ( lat[i] ) * ( lat[Neigh[j]] ); }
    
    // Field term: interaction with external field h
    deltaE += 2. * h * lat[i];
    
    if (deltaE < 0) { w = 1.0; }
    else            { w = exp(-b*deltaE); }
    
    if (drand01() < w ) { lat[i] = -lat[i]; }
}

void simulation::run(int steps, double b) {
    for(int i=0;i<steps;i++) {
    int site = (int) N*drand01();
    MCstep(site,b);
    }
}

void simulation::tempsweep(int steps, double bs, double bf, int bres) {
    std::cout<<"b\t"<<"m\t"<<"m2\t"<<"e\t"<<"e2\t"<<"RR\n";;

    for(int i=0;i<bres;i++) {
    double b = bs + (bf-bs)*i/bres;
    run(steps,b);
    std::cout<<b<<"\t"<<mag()<<"\t"<<mag2()<<"\t"<<ener()<<"\t"<<ener2()<<"\n";

    }

}

void simulation::writeLat() {

    for(int i=0;i<L;i++) {
        for(int j=0;j<L;j++) {
            if ( lat[i*L+j] == 1 ) 
                std::cout<<"1";
            else
                std::cout<<"0";
        }
        std::cout<<"\n";
    }
}

/*
* Write lattice to shared memory buffer
* Format: [8 bytes step_count][4 bytes L][N bytes lattice as int8]
*/
void simulation::writeToShm(char* shm_ptr, long step_count) {
    // Write header: step count (8 bytes) + L (4 bytes)
    memcpy(shm_ptr, &step_count, sizeof(long));
    memcpy(shm_ptr + sizeof(long), &L, sizeof(int));
    
    // Write lattice as int8: +1 or -1
    char* lat_ptr = shm_ptr + sizeof(long) + sizeof(int);
    for(int i = 0; i < N; i++) {
        lat_ptr[i] = (char)lat[i];  // lat[i] is +1 or -1
    }
}

// Global for signal handler cleanup
static int shm_fd_global = -1;
static const char* shm_name_global = "/ising_lattice";

void cleanup_handler(int sig) {
    if(shm_fd_global >= 0) {
        close(shm_fd_global);
        shm_unlink(shm_name_global);
    }
    exit(0);
}
 
int main(int argc, char* argv[]) {
        
    int L = atoi( argv[1] );                              // lattice size
    double beta = atof( argv[2] );                        // inverse temperature
    double J = argc > 3 ? atof(argv[3]) : 1.0;           // coupling (default 1)
    double h = argc > 4 ? atof(argv[4]) : 0.0;           // field (default 0)
    int N = L*L;    // number of sites
    int N2 = N*N;   // square of number of sites
    int NN = 4;     // number of nearest neighbors

    simulation Sim(L, J, h);
    Sim.init();
    
    // Setup shared memory
    const char* shm_name = "/ising_lattice";
    size_t shm_size = sizeof(long) + sizeof(int) + N;  // step + L + lattice
    
    int shm_fd = shm_open(shm_name, O_CREAT | O_RDWR, 0666);
    if(shm_fd < 0) {
        perror("shm_open failed");
        return 1;
    }
    
    if(ftruncate(shm_fd, shm_size) < 0) {
        perror("ftruncate failed");
        return 1;
    }
    
    char* shm_ptr = (char*)mmap(NULL, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if(shm_ptr == MAP_FAILED) {
        perror("mmap failed");
        return 1;
    }
    
    // Setup cleanup on exit
    shm_fd_global = shm_fd;
    signal(SIGINT, cleanup_handler);
    signal(SIGTERM, cleanup_handler);
    
    std::cout << "Ising simulation: L=" << L << ", β=" << beta << ", J=" << J << ", h=" << h << std::endl;
    std::cout << "Shared memory: " << shm_name << " (" << shm_size << " bytes)" << std::endl;
    
    long step_count = 0;
    int steps_per_update = N;  // N steps = 1 sweep on average

    // Run indefinitely until user stops
    while(true) {
        Sim.run(steps_per_update, beta);
        step_count += steps_per_update;
        
        // Write to shared memory
        Sim.writeToShm(shm_ptr, step_count);
    }
    
    // Cleanup (unreachable, but good practice)
    munmap(shm_ptr, shm_size);
    close(shm_fd);
    shm_unlink(shm_name);
    
    return 0;
}
