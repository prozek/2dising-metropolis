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
    simulation(int inL);
    ~simulation();

    double mag();
    double mag2();
    double ener();
    double ener2();
    
    int* genNeigh(int i);
    
    void check();
    void init();
    void MCstep(int site, double b);
    void run(int steps, double b);
    void tempsweep(int steps, double bs, double bf);

private:
    int L;
    int N;
    int NN;
    int* lat;
};

simulation::simulation() {}

simulation::simulation(int inL) {
    L = inL;
    N = L*L;
    NN = 4;
    lat = new int[N];
}

simulation::~simulation() {}

/*
* Generate positions of nearest neighbors on given site, using time-expensive method.
*/
int* simulation::genNeigh(int i) {
    
    int L = 4;
    int N = L*L;
    int NN = 4;

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
	for(int i=0;i<N;i++) {			
		for(int j=0;j<NN;j++) {
		res += (double) lat[i]*lat[genNeigh(i)[j]]; 
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


void simulation::check() {
    int t = 5;
    for(int i=0;i<4;i++)
        std::cout<<genNeigh(t)[i]<<"\n";
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
    
    for (int j=0;j<NN;j++) { deltaE += 2. * ( lat[i] ) * ( lat[genNeigh(i)[j]] ); }
    
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

void simulation::tempsweep(int steps, double bs, double bf) {

}
 
int main(int argc, char* argv[]) {
        
    int L = 4; //atoi( argv[1] );     // lattice size
    int N = L*L;    // number of sites
    int NN = 4;     // number of nearest neighbors

    simulation Sim(L);

    return 0;
}