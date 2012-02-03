#include "cepmods.h"

void hjmap_mod(int nfiles, float* tdata, float* vdata, float* hdata, float* jdata, int freq, int delay, int nodect, MapStruct* hjs)
{
//    *hjs = (float*)malloc(nodect*sizeof(float));
//    if(*hjs == NULL)
//	all_abort("Failed to allocate memory for hjmap_mod.");

    for(int node = 0; node < nodect; node++) {
        hj_mod(tdata, &vdata[node*nfiles], &hdata[node*nfiles], &jdata[node*nfiles], nfiles, freq, delay, &hjs[node]);
    }
}
