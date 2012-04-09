#include "cepmods.h"
#include <math.h>

// Computes time-domain metrics used by Seth Weinberg for HFAC pulse analysis
// Written by Kelly Chang
void mv_mod(float* tdata, float* vdata, int ndata, int freq, MapStruct* vmaxs, MapStruct* vmins, MapStruct* velevs, MapStruct* deltavs)
{
    if(ndata <= 0){
	all_abort("No data provided to mv_mod.");
	exit(1);
    }

#ifdef DEBUG
    cerr << "Received count of ndata: " << ndata << endl;
#endif

    float timeframe = tdata[ndata-1] - tdata[0];
    float sampfreq = (ndata-1)/timeframe;
    int steps = sampfreq/freq; // number of time steps per cycle

    // Find minimum Vm
    float minv = vdata[0];
    int imin = 0;
    for(int i = 1; i < ndata; i++){
        if(vdata[i] < minv){
            minv = vdata[i];
            imin = i;
        }
    }

    // Find maximum Vm and Velev before minimum Vm
    float maxv = minv;
    float velev = vdata[imin];
    int istart = imin-4*steps+1;
    int istop = imin;
    if(istart < 0){
        istart = 0;
        istop = 4*steps;
        velev = 0;
    }
    for(int i = istart; i < istop; i++){
        velev += vdata[i];
        if(vdata[i] > maxv)
            maxv = vdata[i];
    }
    velev /= float(4*steps);

    // Compute delta Vm, Vmin, and Vmax
    float deltav = maxv - minv;
    float vmin = velev - deltav/2;
    float vmax = velev + deltav/2;

#ifdef DEBUG
    cerr << "Vmax = " << vmax << ", Vmin = " << vmin << ", Velev = " << velev << ", and deltaV = " << deltav << endl;
#endif

    vmaxs->mapval = vmax;
    vmins->mapval = vmin;
    velevs->mapval = velev;
    deltavs->mapval = deltav;
}
