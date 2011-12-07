#include "cepmods.h"
#include <fftw3.h>
#include <math.h>

// Requires fftw3
// based on code by Hermenegild Arevalo in the lab of Natalia Tryanova
// modified by Kelly Chang to compute frequency domain refractoriness
void fdr_mod(float* tdata, float* vdata, int ndata, int freq, MapStruct* cps, MapStruct* dcps, MapStruct* hfps, MapStruct* tps)
{
    if(ndata <= 0){
	all_abort("No data provided to fdr_mod.");
	exit(1);
    }

#ifdef DEBUG
    cerr << "Received count of ndata: " << ndata << endl;
#endif

    fftw_complex *in, *out;
    fftw_plan p;
    float timeframe = tdata[ndata-1] - tdata[0];
    float sampfreq = (ndata-1)/timeframe;
    float f_conduction[2] = {1, 30}; // range of conduction frequencies (Hz)
    float f_hf[2] = {freq-5, freq+5}; // range of around the HFAC frequency (Hz)

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
    double cond_pow = 0;
    double dc_pow = 0;
    double hf_pow = 0;
    double tot_pow = 0;
    
    for(int j = 0; j < (ndata/2); j++){
	val = (out[j][0]*out[j][0]+out[j][1]*out[j][1])/ndata; //power spectrum
        tot_pow += val;
        if(j*sampfreq/ndata < f_conduction[0])
            dc_pow += val;
        else if(j*sampfreq/ndata >= f_conduction[0] && j*sampfreq/ndata <= f_conduction[1])
            cond_pow += val;
        else if(j*sampfreq/ndata >= f_hf[0] && j*sampfreq/ndata <= f_hf[1])
            hf_pow += val;
    }

#ifdef DEBUG
    cerr << "Total power = " << tot_pow << ", Conduction power = " << cond_pow << ", DC power = " << dc_pow << ", HF power = " << hf_pow << ", Samp freq = " << sampfreq << endl;
#endif
    
    cps->mapval = cond_pow; // power in the conduction range
    dcps->mapval = dc_pow; // power in the DC range (frequencies less than conduction)
    hfps->mapval = hf_pow; // power around the HFAC frequency
    tps->mapval = tot_pow; // total power

    free(in);
    free(out);
}
