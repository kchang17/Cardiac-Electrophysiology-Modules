#include "cepmods.h"
#include <math.h>

// 
void mv_mod(float* tdata, float* vdata, float* hdata, float* jdata, int ndata, int freq, int delay, MapStruct* mvs)
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
    int steps = sampfreq/freq;
    int nave = (ndata-1)/steps+1;

    float* hjdata = new float [ndata];
    // Make hj vector
    for(int i = 0; i < ndata; i++){
        hjdata[i] = hdata[i] * jdata[i];
    }

    float* tave = new float [nave];
    float* vave = new float [nave];

    // Get time-averaged Vm trace over each cycle
    for(int i = 0; i < nave; i++){
        tave[i] = tdata[(i+1)*steps];
        for(int j = 0; j < steps; j++){
            vave[i] += vdata[i*steps+j+1];
        }
        vave[i] /= steps;
    }

    // Find the maximum derivative of the averaged trace
    float* dvdt = new float [nave-1];
    float maxdvdt = dvdt[0];
    for(int i = 0; i < nave-1; i++){
        dvdt[i] = (vave[i+1]-vave[i])/(tave[i+1]-tave[i]);
        if(dvdt[i] > maxdvdt){
            maxdvdt = dvdt[i];
        }
    }

    // Find activations
    bool* hasact = new bool [nave-1];
    for(int i = 0; i < nave-1; i++){
        if(maxdvdt > 0 && dvdt[i] >= 0.2*maxdvdt){
            hasact[i] = true;
        }
        else if(maxdvdt <=0 && dvdt[i] <= 0.2*maxdvdt){
            hasact[i] = true;
        }
    }

    // Eliminate adjacent activation detections
    for(int i = 0; i < nave-2; i++){
        if(hasact[i] && hasact[i+1]){
            hasact[i] = false;
        }
    }

    // Eliminate activations that occur within the delay
    int buff = delay*freq/1000;
    int i = 0;
    while(i < nave-1){
        if(hasact[i]){
            for(int j = 0; j < buff; j++){
                i++;
                if(i == nave-1){
                    break;
                }
                hasact[i] = false;
            }
        }
        i++;
    }

    // Find maximum vm one cycle before maximum h*j before each activation
    float tstart = tdata[0];
    float tstop;
    float mvmean = 0;
    int ctr = 0;
    for(int i = 0; i < nave-1; i++){
        if(hasact[i]){
            tstop = tave[i];
//            if(tstart>tdata[0] || tstop >= tstart+float(delay)*0.001){
            int istart = (tstart-tdata[0])*sampfreq+1;
            int istop = (tstop-tdata[0])*sampfreq;
            int ihjmax = istart;
            for(int j = istart+1; j < istop; j++){
                if(hjdata[j]>hjdata[ihjmax]){
                    ihjmax = j;
                }
            }
            float vmmax = vdata[ihjmax];
            for(int j = ihjmax-steps; j < ihjmax; j++){
                if(j>=0 && vdata[j]>vmmax){
                    vmmax = vdata[j];
                }
            }
            mvmean += vmmax;
            ctr++;
//            }
            tstart = tstop;
        }
    }
    if(ctr==0){
        cerr << "No activations detected." << endl;
    }
    mvmean/=float(ctr);

#ifdef DEBUG
    cerr << "mv metric = " << mvmean << ", computed over " << ctr << " activations." << endl;
#endif

    mvs->mapval = mvmean; // mv metric
}
