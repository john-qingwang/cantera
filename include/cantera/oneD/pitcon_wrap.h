#ifndef pitcon_wrap_h
#define pitcon_wrap_h

extern "C" {

void c_f(int*,double*,int*,double*,double*);
void c_pitcon(void (*c_f)(int*,double*,int*,double*,double*),double*,int*,int*,
              int*,int*,int*,double*,int*,double*,int*,double*,double*);
}

#endif
