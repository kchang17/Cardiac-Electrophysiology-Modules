#include "cepmods.h"

void mvmap_mod(int nfiles, float* tdata, float* vdata, float* hdata, float* jdata, int freq, int delay, int nodect, MapStruct* mvs)
{
//    *mvs = (float*)malloc(nodect*sizeof(float));
//    if(*mvs == NULL)
//	all_abort("Failed to allocate memory for mvmap_mod.");

    for(int node = 0; node < nodect; node++) {
        mv_mod(tdata, &vdata[node*nfiles], &hdata[node*nfiles], &jdata[node*nfiles], nfiles, freq, delay, &mvs[node]);
    }
}
