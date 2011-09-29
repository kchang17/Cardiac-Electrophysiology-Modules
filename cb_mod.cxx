#include "cepmods.h"
#include <fftw3.h>
#include <math.h>

// Requires fftw3
// based on code by Hermenegild Arevalo in the lab of Natalia Tryanova
// modified by Kelly Chang to compute metric for conduction block
void cb_mod(float* tdata, float* vdata, int ndata, MapStruct* cbs)
{
    if(ndata <= 0){
	all_abort("No data provided to cb_mod.");
	exit(1);
    }

#ifdef DEBUG
    cerr << "Received count of ndata: " << ndata << endl;
#endif

    fftw_complex *in, *out;
    fftw_plan p;
    float timeframe = tdata[ndata-1] - tdata[0];
    float sampfreq = (ndata-1)/timeframe;
    double f_conduction[2] = {1, 30}; // range of conduction frequencies (Hz)

    in = (double (*)[2])fftw_malloc(sizeof(fftw_complex) * ndata);
    out = (double (*)[2])fftw_malloc(sizeof(fftw_complex) * ndata);
    if(in == NULL || out == NULL)
	all_abort("Failed to allocate arrays needed for FFT.");

    p = fftw_plan_dft_1d(ndata, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for(int tstep = 0; tstep < ndata; tstep++){
	in[tstep][0] = (double)vdata[tstep];
	in[tstep][1] = 0;
    }
    
    fftw_execute(p); // perform fft
    
    double val = 0;
    double tot_area = 0;
    double cond_area = 0;
    double norm_area = 0;
    
    for(int j = 0; j < (ndata/2); j++){
	val = sqrt(out[j][0]*out[j][0]+out[j][1]*out[j][1]); //frequency spectrum
        tot_area += val;
        if(j*sampfreq/ndata >= f_conduction[0] && j*sampfreq/ndata <= f_conduction[1])
            cond_area += val;
    }

    norm_area = cond_area/tot_area; // normalized area in freq range of conduction

#ifdef DEBUG
    cerr << "Total area = " << tot_area << ", Conduction area = " << cond_area << ", Samp freq = " << sampfreq << endl;
#endif
    
    cbs->mapval = 100*(1-norm_area); // percent area outside freq range of conduction

    free(in);
    free(out);
}
