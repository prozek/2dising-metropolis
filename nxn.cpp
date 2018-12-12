#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>


class simulation {
public:
    simulation();
    simulation(int inL);
    //~simulation();


    int L;
    /*int N = L*L;
    int NN = 4;
    int* lat = new int[N];
*/
};

simulation::simulation() {
    std::cout<<"Constructor works!";
}

simulation::simulation(int inL) {
    L = inL;
    std::cout<<"Constructor works!";
}

/*
* Generate positions of nearest neighbors on given site.
* char spt when assigned value 'M' generated positions are assigned to the table Neigh. When value 'T' is assigned neigbor list is generated each time function is called.
*
*/
int* genNeigh(char spt, int x, int y) {
    
    int L = 4;
    int N = L*L;
    int NN = 4;

        int* Neigh = new int[NN];
        int i = x + L*y;
        
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


double mag(int *lat, int NN, int N) {
	double res=0;
	for(int i=0;i<N;i++) {
		res += (double) lat[i]; // *(lat+i);	
    }
    return res ;// (double) N;
}

double mag2(int *lat, int NN, int N) {
	double res=0;
	for(int i=0;i<N;i++) {
	for(int j=0;j<N;j++) {
		res += (double) lat[i]*lat[j]; // ( (*(lat+i))*(*(lat+j)) );	
	}}	
	return res ;// (double) (N*N);
}

double energy2(int **Neigh, int *lat, int NN, int N) {
	double res=0;	
	for(int i1=0;i1<N;i1++) {
	for(int j1=0;j1<N;j1++) {
		for(int k1=0;k1<NN;k1++) { 
		for(int l1=0;l1<NN;l1++) {	
			res += (double) lat[i1]*lat[Neigh[i1][k1]] * lat[j1]*lat[Neigh[j1][l1]];
			}}
		}}
	return res;
}


double energy(int **Neigh, int *lat, int NN, int N) {
	double res=0;
	for(int i=0;i<N;i++) {			
		for(int j=0;j<NN;j++) {
		res += (double) lat[i]*lat[Neigh[i][j]]; 
		}
	}
	return res; // (double) N; 
}

double drand01() {
    return double( rand() )/ double( RAND_MAX );
}
 
int main(int argc, char* argv[]) {
        
    int L = 4; //atoi( argv[1] );     // lattice size
    int N = L*L;    // number of sites
    int NN = 4;     // number of nearest neighbors

    simulation Sim = simulation(4);
    Sim.L = 4;
    std::cout<<Sim.L;

/*  int mcsmax = 1000;

    int lattice[N];
    int** Neigh = new int*[N];
    for(int i =0;i<N;++i) {
        Neigh[i] = new int[NN]; }

    int t, tres=100;
    double beta, endt=4, begt= 0.1;
    int cutstep=10;

    int* pos = genNeigh('a',1,1);
    std::cout<<"\n";
    for(int i=0;i<4;i++)
        std::cout<<pos[i]<<"\n";

    for (int i = 0; i < N; i++) {
        double rd = drand();
        if (rd < 0.5) { lattice[i] = 1; } else { lattice[i] = -1; }
    }

	for(t=0; t <= tres; t++) {      // temperature loop
	double mtot=0;
	double mfluc=0;
	double m2=0;
	double enrg=0;
	double enrg2=0;


for (int mcs = 1; mcs <= 1000; mcs++) {   // thermalization
        double deltaE;
            
        for (int i = 0; i < N; i++) {	
      		deltaE=0.0;
 	        double w;
            for (int j=0; j<NN; j++) { deltaE += 2.* ( lattice[i] ) * ( lattice[Neigh[i][j]] ); }
            
	        double rd = drand();
	        beta = 1. / (double) ( begt + t * (endt-begt) / (double) (tres) );

	        if (deltaE < 0) { w = 1.0 ; }
	        else            { w = exp(-beta*deltaE) ; }

      	    if (rd<w) { lattice[i] = -lattice[i] ; }
            }
    }

    for (int mcs = 1; mcs <= mcsmax; mcs++) {   // MC loop
        double deltaE;
            
        for (int i = 0; i < N; i++) {	
      		deltaE=0.0;
 	        double w;
            for (int j=0; j<NN; j++) { deltaE += 2.* ( lattice[i] ) * ( lattice[Neigh[i][j]] ); }
            
	        double rd = drand();
	        beta = 1. / (double) ( begt + t * (endt-begt) / (double) (tres) );

	        if (deltaE < 0) { w = 1.0 ; }
	        else            { w = exp(-beta*deltaE) ; }

      	    if (rd<w) { lattice[i] = -lattice[i] ; }
            }

	if( mcs%cutstep == 0 ) {
		mtot+=	cutstep*fabs( mag(lattice,NN,N)  );
		m2+=	cutstep*fabs( mag2(lattice,NN,N) );
		enrg+=	cutstep*( energy(Neigh,lattice,NN,N)   );
		enrg2+=	cutstep*( energy2(Neigh,lattice,NN,N)  );
		}
            } //mc steps 
	 
        //std::cout<<2/beta<<" , "<<mtot/mcsmax<<" , "<<(m2-mtot*mtot/mcsmax)/mcsmax<<" , "<<enrg/mcsmax<<" , "<<(enrg2-enrg*enrg/mcsmax)/mcsmax<<"\n";
        
        } //b steps	
*/
//    writeLatt(lattice,N,L);
//    std::cout<<"\n";
    
//    for(int i=0;i<N;i++) { std::cout<<i<<" ; "; 
//    for(int j=0;j<NN;j++) { std::cout<<Neigh[i][j]<<" , "; }
//    std::cout<<"\n";}
//    writevec(lattice,N,L);
    return 0;
}
