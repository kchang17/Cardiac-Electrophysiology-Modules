#include "cepmods.h"

void cbmap_mod(int nfiles, float* tdata, float* vdata, int nodect, MapStruct* cbs)
{
//    *cbs = (float*)malloc(nodect*sizeof(float));
//    if(*cbs == NULL)
//	all_abort("Failed to allocate memory for cbmap_mod.");

    for(int node = 0; node < nodect; node++) {
        cb_mod(tdata, &vdata[node*nfiles], nfiles, &cbs[node]);
    }
}
