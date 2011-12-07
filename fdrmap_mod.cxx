#include "cepmods.h"

void fdrmap_mod(int nfiles, float* tdata, float* vdata, int freq, int nodect, MapStruct* cps, MapStruct* dcps, MapStruct* hfps, MapStruct* tps)
{
    for(int node = 0; node < nodect; node++) {
        fdr_mod(tdata, &vdata[node*nfiles], nfiles, freq, &cps[node], &dcps[node], &hfps[node], &tps[node]);
    }
}
