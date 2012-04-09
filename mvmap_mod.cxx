#include "cepmods.h"

void mvmap_mod(int nfiles, float* tdata, float* vdata, int freq, int nodect, MapStruct* vmaxs, MapStruct* vmins, MapStruct* velevs, MapStruct* deltavs)
{
//    *vmaxs = (float*)malloc(nodect*sizeof(float));
//    *vmins = (float*)malloc(nodect*sizeof(float));
//    *velevs = (float*)malloc(nodect*sizeof(float));
//    *deltavs = (float*)malloc(nodect*sizeof(float));
//    if(*vmaxs == NULL || *vmins == NULL || *velevs == NULL || *deltavs == NULL)
//	all_abort("Failed to allocate memory for mvmap_mod.");

    for(int node = 0; node < nodect; node++) {
        mv_mod(tdata, &vdata[node*nfiles], nfiles, freq, &vmaxs[node], &vmins[node], &velevs[node], &deltavs[node]);
    }
}
