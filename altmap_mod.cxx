#include "cepmods.h"

void altmap_mod(int nfiles, float* tdata, float* vdata, int nodect, MapStruct* apds, MapStruct* alts, int bcl, int t, int dim_t)
{

#ifdef DEBUG
    cerr << "Mapping " << nodect << " nodes from " << nfiles << " files."<< endl;
#endif

    for(int node = 0; node < nodect; node++){
#ifdef DEBUG
//	cerr << "Passing node " << node << " to alt_mod(). " << endl;
#endif
	alt_mod(tdata, &vdata[node*nfiles], nfiles, &apds[node], &alts[node], bcl, t, dim_t);
#ifdef DEBUG
//	cerr << (*apds+node)->time << endl;
#endif
    }
}
