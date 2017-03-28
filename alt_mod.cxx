/* 
 * Kelly Chang
 * Alternans/APD90 Module
 * 
 *  alt_mod.cxx - finds an alternans and APD90 in a time/voltage array
 *
 */

#include "cepmods.h"
#include <cmath>

/* data is a 2D array arranged in columns */
/* column 1 is the time, column 2 is vm */
void alt_mod(float* tdata, float* vdata, int ndata, MapStruct* apddata, MapStruct* altdata, int bcl, int t, int dim_t)
{
    
#ifdef DEBUG
//    cerr << "Kelly's APD/alternans data: "<< endl;
#endif
    float dvdt;
    float d2vdt2;
    float maxv;
    float maxdvdt;
    float maxd2vdt2;
    int status;
    int maxvtime;
    int i;
    int j;
    int nbeats;
    int bcl_ndata;
    bool block;
    int diptime;

    if (bcl <= 100)
        cerr << "BCL should be above 100 ms." << endl;
    nbeats = (ndata-1)*(dim_t+50)/(t-1)/bcl-2;
    if (nbeats < 1)
        cerr << "Need at least 3 beats to get DI-APD pair." << endl;
    if (fmod(float(bcl)*float(t-1),float(dim_t)) != 0)
        cerr << "BCL is not a multiple of fac_t=(t-1)/dim_t!" << endl;
    bcl_ndata = bcl*(t-1)/dim_t;

    int acttime[nbeats+2];
    int notchtime[nbeats+2];
    int pII_time[nbeats+2];
    int pIV_time[nbeats+2];
    float pIV_amp[nbeats+2];
    float pII_amp[nbeats+2];
    for (i = 0; i < nbeats+2; i++){
        maxdvdt = 0;
        acttime[i] = -1;
        status = 0;
        pIV_amp[i] = vdata[i*bcl_ndata];
        maxv = -100;
        maxvtime = -1;
        maxd2vdt2 = -1000;
        block = false;
        diptime = -100;

        int last_pt;
        if( i == nbeats+1 )
            last_pt = ndata;
        else
            last_pt = (i+1)*bcl_ndata;
        for (j = i*bcl_ndata+1; j < last_pt; j++){
            // Get phase IV amplitude
            if (status == 0 && pIV_amp[i] > vdata[j]){
                pIV_amp[i] = vdata[j];
                pIV_time[i] = j;
            }
            // Get activation time
            dvdt = (vdata[j] - vdata[j-1])/(tdata[j] - tdata[j-1]);
            if (dvdt > maxdvdt){
                maxdvdt = dvdt;
                acttime[i] = j-1;
                maxv = vdata[j];
                maxvtime = j;
                status = 1;
                block = true;
            }else if(status == 1 && dvdt < 0){
                // if there's a blip before threshold crossing, reset max dvdt
                maxdvdt = 0;
                acttime[i] = -1;
                maxv = -100;
                maxvtime = -1;
            }
            // Make sure Vm crosses -30 mV
            if (block && vdata[j] >= -0.030) {
                block = false;
                status = -1;
            }
            // Find peak
            if (status == -1){
                if (vdata[j] > maxv){
                    maxv = vdata[j];
                    maxvtime = j;
                }
                else if (dvdt < 0){
                    notchtime[i] = maxvtime;
                    pII_amp[i] = maxv;
                    pII_time[i] = maxvtime;
                    status = -2;
                }
            }
            // Find notch
            if (!block && status == -2){
                if (maxvtime - acttime[i] > 5){ // this was 4 in MATLAB code?
                    // If peak is delayed, let phase II be the peak
                    notchtime[i] = maxvtime-1;
                    pII_amp[i] = maxv;
                    pII_time[i] = maxvtime;
                    status = -3;
                }
                else if (j <= maxvtime+100 && dvdt >= 0){
                    // Check if there's an actual notch within 100 ms
                    notchtime[i] = j-1;
                    pII_amp[i] = vdata[j];
                    pII_time[i] = j;
                    diptime = j-1;
                    status = -3;
                }
                else if (diptime == -100 && j > maxvtime+1 && j < maxvtime+100){
                    d2vdt2 = (vdata[j+1]-vdata[j])/(tdata[j+1]-tdata[j])-(vdata[j]-vdata[j-1])/(tdata[j]-tdata[j-1]);
                    if (d2vdt2 > maxd2vdt2){
                        maxd2vdt2 = d2vdt2;
                        notchtime[i] = j;
                        pII_amp[i] = vdata[j];
                        pII_time[i] = j;
                    }
                }
            }
            // Check for extra blips (notch should be the last dip)
            if (!block && status == -3 && j<diptime+5 && dvdt < 0){
                diptime = j;
                status = -2;
            }
            // Find phase II amplitude
            if (!block && (status == -2 || status == -3)){
                if (notchtime[i] - maxvtime > 50){
                    // If notch is delayed, set it to be right after the peak
                    notchtime[i] = maxvtime+1;
                    pII_amp[i] = vdata[maxvtime+1];
                    pII_time[i] = maxvtime+1;
                    status = -4; // should this be -3 to match MATLAB code?
                }
                if (j <= notchtime[i]+100){
                    // Check for maximum within 100 ms of notch
                    if (vdata[j] > pII_amp[i]){
                        pII_amp[i] = vdata[j];
                        pII_time[i] = j;
                    }
                }
                else{
                    status = -4;
                }
            }
        }
        if (block){
            // no activation detected
            apddata->mapval = -1;
            altdata->mapval = -1;
            return;
        }
    }

    // Calculate APD90
    float apdsum;
    float altsum;
    float apd90;
    float apd90_prev;
    float pIV_mean;
    float v90;
    int apd90_time;
    int ctr;
    int factor;
    apdsum = 0;
    altsum = 0;
    ctr = 0;
    factor = 1;
    for (i = 0; i < nbeats+1; i++){
        // Get phase IV amplitude
        pIV_mean = (pIV_amp[i]+pIV_amp[i+1])/2;
        v90 = pIV_mean+0.1*(pII_amp[i]-pIV_amp[i]);
        // Find APD90
        apd90_time = acttime[i+1];
        for (j = pII_time[i]; j < acttime[i+1]; j++){
            if (j < apd90_time && vdata[j] < v90){
                apd90_time = j;
                apd90 = tdata[apd90_time]-tdata[acttime[i]];
                if (i > 0){
                    apdsum += apd90;
                    altsum += factor*(apd90-apd90_prev);
                    ctr++;
                    factor *= -1;
                }
                apd90_prev = apd90;
            }
            if (j < apd90_time){
                // Check for funny stuff
                dvdt = (vdata[j+1]-vdata[j])/(tdata[j+1]-tdata[j]);
                if (dvdt > 0.6){
                    apddata->mapval = -1;
                    altdata->mapval = -1;
                    return;
                }
            }
            else if (j < pIV_amp[i] && vdata[j] > -40){
                // Check for activations after first AP
                apddata->mapval = -1;
                altdata->mapval = -1;
                return;
            }
        }
    }
    // Return average APD90 and average alternans magnitude (with sign)
    apddata->mapval = apdsum/ctr;
    altdata->mapval = altsum/ctr;
}

