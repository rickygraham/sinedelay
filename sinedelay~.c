#include "m_pd.h"
#include "math.h"
#include "stdlib.h"
#define TWOPI 6.2831853072
#define XTRASAMPS 4

/* implement later - this external assumes block of 64 samples for now */
extern int ugen_getsortno(void); //in d_ugen.c returns dsp sort number 

static t_class *sinedelay_class;

typedef struct _sinedelay
{
  t_object obj;
  t_sample x_f;
  long table_length;
  float *wavetable; 
  float *window; 
  float sr, phase, si, si_w, phase_w;
  long wavetable_bytes;
  long window_bytes;
  int c;

 /*----for delay line----*/
  long delay_length; // length of the delay line in samples
  long delay_bytes; // length of delay line in bytes
  float maximum_delay_time, feedback, *delay_line;

  long write_index; // write point in delay line

  /*---for crossfade-----*/
  int dirty;
  float *wavetable_old; 
  float s_fraction;
  float xfade_samples;
  float xfade_countdown;
  float fdelay;


  /*--------4-point interpolation-----*/
  float zerodel;
  int c_phase;
  int c_n;
  float *c_vec;

} t_sinedelay;


float biclip(float n)
{
    return (n < -1) ? -1.0 : (n > 1) ? 1.0 : n;
}


float uniclip(float n, float o) //determine clip point
{
    return (n < 0) ? 0 : (n > o) ? o : n;
}


float softclip(float n)
{
    return (n < -0.5) ? tanh(n) : (n > 0.5) ? tanh(n) : n;
}



int modulo(int a, int b)
{
	//n = input; a = max value; 
    return (a/b)*b + a%b;
}


float LinearInterpolate(float y1, float y2, float mu) {
	return(y1*(1-mu)+y2*mu);
	}


void sinebasic(t_sinedelay *x) {
	
	if(x->xfade_countdown <= 0)x->xfade_countdown=0;

	for (int i = 0; i < x->table_length; i++) {
	x->wavetable_old[i] = x->wavetable[i];
	} 

	x->dirty=1;
	x->xfade_countdown = x->xfade_samples; 

	for(int i = 0; i < x->table_length; i++) {
	x->wavetable[i] = sin(TWOPI * (float)i / (float)x->table_length);
	} 


	x->dirty=0;

	/* window for am - modify per waveform*/
	float tp; 
	tp = (float)x->table_length * 2;
	for(int a = 0; a < tp; a++) {
	x->window[a] = sin(TWOPI * (float)a / tp);
	}    

}


void squarewave(t_sinedelay *x) {

	if(x->xfade_countdown <= 0)x->xfade_countdown=0;

	for (int i = 0; i < x->table_length; i++) { 
	x->wavetable_old[i] = x->wavetable[i];
	} 

	x->dirty=1;

	x->xfade_countdown = x->xfade_samples;

	for(int i = 0; i < x->table_length; i++) {
    float sqr = 0;
    if(i > x->table_length/2){
    sqr = 1.0;
    } else {
    sqr = -1.0;
    }
    x->wavetable[i] = sqr;
	}

	x->dirty=0;

  /* window for am - modify per waveform */
  float tp; 
	tp = (float)x->table_length * 2;
	for(int a = 0; a < tp; a++) {
	x->window[a] = sin(TWOPI * (float)a / tp);
	}
  
}

void *sinedelay_new(void)
{

  t_sinedelay *x = (t_sinedelay *) pd_new(sinedelay_class);
  
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //delaytime
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //feedback
  inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal")); //am
  outlet_new(&x->obj, gensym("signal")); //dry out
  outlet_new(&x->obj, gensym("signal")); //wet out

  float init_freq = 440.0;
  x->table_length = 65536;
  
  x->wavetable_bytes = x->table_length * sizeof(float);
  x->wavetable = (float *) getbytes(x->wavetable_bytes);
  x->wavetable_old = (float *) getbytes(x->wavetable_bytes);
  
  x->window_bytes = (x->table_length * 2) * sizeof(float);
  x->window = (float *) getbytes(x->window_bytes);
  
  x->phase = 0; x->phase_w = 0; x->c = 0;
  x->sr = sys_getsr(); 
  x->si = init_freq * ((float)x->table_length/x->sr);

   /* delay line */
  x->feedback = 0.5;
  x->fdelay = 150 * x->sr / 1000;
  x->zerodel = 64; //4-point interpolation - v1. assumes block is 64 samples
  x->maximum_delay_time = 10;
  x->delay_length = x->sr * x->maximum_delay_time;
  x->delay_bytes = x->delay_length * sizeof(float);
  x->delay_line = (float *)getbytes(x->delay_bytes);

  /*-----4 point Interpolation------*/
  x->c_n = 0;
  x->c_vec = getbytes(XTRASAMPS * sizeof(t_sample)); //bytes for 4 sample interpolation 
  
  /*---xfade----*/
  x->xfade_countdown = 1;
  x->xfade_samples = 150 * (x->sr / 1000.0);

  /* clear the delayline */
  for(int i = 0; i < x->delay_length; i++){
    x->delay_line[i] = 0.0;
  }

  x->write_index = 0;

  sinebasic(x);
  return x; //return pointer to object
}


t_int *sinedelay_perform(t_int *w) {
  t_sinedelay *x = (t_sinedelay *) (w[1]);
  float *frequency = (float *)(w[2]);
  float *delaytime = (float *)(w[3]);
  float *feedback = (float *)(w[4]);
  float *am = (float *)(w[5]);
  float *out1 = (float *)(w[6]);
  float *out2 = (float *)(w[7]);
  int n = w[8];
  long iphase, iphase_w;
  float si = x->si;
  float si_w = x->si_w; //for amplitude modulation
  float phase = x->phase;
  float phase_w = x->phase_w; //for amplitude modulation
  long table_length = x->table_length;
  float *wavetable = x->wavetable;
  float *wavetable_old = x->wavetable_old;
  float *window = x->window;
  float sr = x->sr;

  float d_fraction; //interp for reading delay line
  float s_fraction; //for waveform xfade
  float xfade_samples = x->xfade_samples;
  float xfade_countdown = x->xfade_countdown;
  
  /* for delay */
  float *delay_line = x->delay_line;
  long  write_index = x->write_index;
  long  delay_length = x->delay_length; 
  float fdelay = x->fdelay;
  long  idelay;
  float srms = sr / 1000.0;
  float fback = x->feedback;
  float c_wave, o_wave, winsample, waveout, delayout;
  float fn = n-1; //val of 63 assuming a block of 64 samples
  //float *vp = x->delay_line;
  float *bp, *wp = &delay_line[write_index];
  float zerodel = x->zerodel;
  float fourout;

while (n--) {
/* waveform code first - including windowing */
si =   *frequency++ * ((float)table_length/sr); //hertz to samples
si_w = *am++ * ((float)table_length/sr); //table_length val is fine because we only want to read half of the sine
iphase = (int)phase; 
iphase_w = (int)phase_w; 

/* update and constrain phase values*/ 
phase += si; 
phase_w += si_w; 

     while(phase >= table_length) {
     phase -= table_length;
     } 
     while(phase < 0) {
     phase += table_length;
     } 
       while(phase_w >= table_length) {
     phase_w -= table_length;
     }  
     while(phase_w < 0) {
     phase_w += table_length;
     } 

/* assign vals to local vars */
c_wave = wavetable[iphase]; 
o_wave = wavetable_old[iphase];
winsample = window[iphase_w];

/* check flag for wavetable cross-fade */
if(x->dirty){ 
  waveout = o_wave;
  } else if(xfade_countdown > 0.0){ //xfade will be negative if updated with fractional val
  s_fraction = xfade_countdown-- / xfade_samples;
  //post("fraction: %f\n, countdown:%f\n, fadesamps: %f\n", s_fraction, xfade_countdown, xfade_samples);
  waveout = wavetable[iphase] + s_fraction * (wavetable_old[iphase] - wavetable[iphase]); 
  } else {
  waveout = wavetable[iphase];
  }


/*--------------------end of waveform xfade code------------------------*/

/* now, delay line code */
t_sample a, b, c, d, cminusb;

fdelay = *delaytime++ *srms - zerodel; //like delsamps - 64 assuming block size of 64

  while(fdelay > delay_length){
      fdelay -= delay_length;
    }

   while(fdelay < 0){
      fdelay += delay_length;
    }


fdelay += fn; //adding 63 samples back  = n - 1 = 64 - 1
fn = fn - 1.0f; //decremented every sample 

  /*
  ///notes
  nsamps = delaytime in samps
  c_n = nsamps in delwrite update sr writes delay time to c_n in x_cspace in vd~
  c_phase is set to a val of 4 samples
  */

idelay = (int)fdelay;
d_fraction = fdelay - idelay; //like frac in vd~

/* bp is read position - write pointer minus delay - don't read where you write */
bp = wp - idelay; 

/* four-point interpolation */
if (bp < delay_line + 4) bp += delay_length; //delay_length is max dl time in samps
        d = bp[-3];
        c = bp[-2];
        b = bp[-1];
        a = bp[0];
        cminusb = c-b;
        fourout = b + d_fraction * (
            cminusb - 0.1666667f * (1.- d_fraction) * (
                (d - a - 3.0f * cminusb) * d_fraction + (d + 2.0f*a - 3.0f*b)
            )
        );

/* feedback param input*/
fback = *feedback++;

/* address denormals in feedback param*/
if(fabs(fback) < .000001) { 
  fback = 0.0;
  }

delay_line[write_index++] = waveout + fourout * fback;

/* constrain write index */
if(write_index >= delay_length) {
  write_index -= delay_length;
  }
    
delayout = fourout;
*out1++ = biclip(waveout * winsample * 2);
*out2++ = softclip(delayout * winsample);


} /* ----------- end of block / while loop n-- -------------- */
     x->xfade_countdown = xfade_countdown;
     x->phase = phase; 
     x->phase_w = phase_w;
     x->write_index = write_index;
    return w + 9;
}


void clear(t_sinedelay *x)
{
  for(int i =0; i < x->table_length; i++)
  {
  x->wavetable[i] = 0.0;
  x->wavetable_old[i] = 0.0;
  }

  for(int i = 0; i < x->delay_length; i++){
    x->delay_line[i] = 0.0;
  }

}


void xfade(t_sinedelay *x, float f)
{
if(f>0.0)x->xfade_samples = f * (x->sr / 1000.0);
}


void sinedelay_dsp(t_sinedelay *x, t_signal **sp)
{
    dsp_add(sinedelay_perform, 8, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
}


void sinedelay_free(t_sinedelay *x)
{
  t_freebytes(x->wavetable, x->wavetable_bytes);
  t_freebytes(x->wavetable_old, x->wavetable_bytes);
  t_freebytes(x->window, x->window_bytes);
  t_freebytes(x->delay_line, x->delay_bytes);
  t_freebytes(x->c_vec, x->c_n + XTRASAMPS * sizeof(t_sample)); //free 4 bytes for 4-point interpolation
}


void sinedelay_tilde_setup(void)
{
  sinedelay_class = class_new(gensym("sinedelay~"), (t_newmethod)sinedelay_new, (t_method)sinedelay_free, sizeof(t_sinedelay), CLASS_DEFAULT, A_GIMME, 0);
  class_addmethod(sinedelay_class, (t_method)sinedelay_dsp, gensym("dsp"),0);
  CLASS_MAINSIGNALIN(sinedelay_class, t_sinedelay, x_f);
  class_addmethod(sinedelay_class, (t_method)sinebasic, gensym("sine"),0);
  class_addmethod(sinedelay_class, (t_method)squarewave, gensym("square"),0);
  class_addmethod(sinedelay_class, (t_method)clear, gensym("clear"),0);
  class_addmethod(sinedelay_class, (t_method)xfade, gensym("xfade"), A_DEFFLOAT, 0);
  post("sinedelay~: show me a sine");
}
