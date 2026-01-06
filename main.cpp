#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>

/*
* Choose a random double from range [0,1]
*/
double drand01() {
    return double( rand() )/ double( RAND_MAX );
}

class simulation {
public:
    simulation();
    simulation(int inL, double coupl);
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
     
private:
    int L;
    int N;
    int NN;
    int* lat;
    double J;
};

simulation::simulation() {}

simulation::simulation(int inL, double coupl) {
    L = inL;
    N = L*L;
    NN = 4;
    lat = new int[N];
    J = coupl;
}

simulation::~simulation() {}

/*
* Generate positions of nearest neighbors on given site, using time-expensive method.
*/
int* simulation::genNeigh(int i) {
    
   int* Neigh = new int[NN];     
        
        if( (i+1)%L == 0 )  Neigh[0]=i+1-L;
        else          Neigh[0]=i+1;

        if( ((i-1)%L == (L-1)) ||  ((i-1)%L == -1) )  Neigh[1]=i-1+L;
        else          Neigh[1]=i-1;
        
        if( (i+L) >= N)  Neigh[2]=(i+L-N)%N;
        else          Neigh[2]=i+L;

        if( (i-L) < 0)  Neigh[3]=(i-L+N)%N;
        else          Neigh[3]=i-L;
        
        return Neigh; 
}

/*
* Optimized inline neighbor lookup - no allocation
*/
void simulation::getNeighbors(int i, int* Neigh) {
    if( (i+1)%L == 0 )  Neigh[0]=i+1-L;
    else          Neigh[0]=i+1;

    if( ((i-1)%L == (L-1)) ||  ((i-1)%L == -1) )  Neigh[1]=i-1+L;
    else          Neigh[1]=i-1;
    
    if( (i+L) >= N)  Neigh[2]=(i+L-N)%N;
    else          Neigh[2]=i+L;

    if( (i-L) < 0)  Neigh[3]=(i-L+N)%N;
    else          Neigh[3]=i-L;
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
    
    for (int j=0;j<NN;j++) { deltaE += 2. * J * ( lat[i] ) * ( lat[Neigh[j]] ); }
    
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
 
int main(int argc, char* argv[]) {
        
    double J = 1; // Heisenberg/Ising coupling term
    int L = atoi( argv[1] );     // lattice size
    double local_temperature = atof( argv[2] );
    int N = L*L;    // number of sites
    int N2 = N*N;   // square of number of sites
    int NN = 4;     // number of nearest neighbors

    simulation Sim(L, J);
    Sim.init();
    
    // Flush output immediately
    std::cout.setf(std::ios::unitbuf);
    
    int chunk_steps = 100;
    int step_count = 0;

    // Run indefinitely until user stops
    for(int i=0; ; i++) {
        Sim.run(chunk_steps, local_temperature);
        step_count += chunk_steps;
        
        // Output frame with step info
        std::cout<<"FRAME "<<step_count<<"\n";
        Sim.writeLat();
        std::cout<<"END\n";
        std::cout.flush();
    }
    
    return 0;
}
