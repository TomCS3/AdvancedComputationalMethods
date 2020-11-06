#include<stdio.h>
/* TIMING CODE BEGIN (We need the following lines to take the timings.) */
#include<stdlib.h>
#include<math.h>
#include <time.h>
clock_t startm, stopm;
#define RUNS 1
#define START if ( (startm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define STOP if ( (stopm = clock()) == -1) {printf("Error calling clock");exit(1);}
#define PRINTTIME printf( "%8.5f seconds used .", (((double) stopm-startm)/CLOCKS_PER_SEC/RUNS));
/* TIMING CODE END */

#define N 10000000

/*function prototypes*/
double f_x(double x);
double pi(long n);

int main(void) {
    /* Declarations */

    /* Code */
    START;               /* Timing measurement starts here */
    /* Code to be written by student, calling functions from here is fine
       if desired
    */
	printf("f(x) = %f\n", f_x(0.1));
	printf("pi(n) = %f\n", pi(N));		

    STOP;                /* Timing measurement stops here */
    PRINTTIME;           /* Print timing results */
    return 0;
}

double f_x(double x){
	return sqrt(1-x*x);
}

double pi(long n){
	double a = -1, b =1, i, h, s, x, PI;
	h = (b-a)/n;
	s = 0.5*f_x(a)+0.5*f_x(b);
	for(i=1; i<=n-1; i++){
		x = a+i*h;
		s = s+f_x(x);
	}
	PI = s*h*2;
	
	return PI;
}	 

