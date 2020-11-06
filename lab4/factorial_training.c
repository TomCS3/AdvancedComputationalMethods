#include <stdio.h>
#include <limits.h>
#include <math.h>


/*function prototypes*/
long maxlong(void);

double upper_bound(long n);

int main(void) {
    long i;

    /* The next line should compile once "maxlong" is defined. */
    printf("maxlong()=%ld\n", maxlong());

    /* The next code block should compile once "upper_bound" is defined. */

    
    for (i=0; i<10; i++) {
        printf("upper_bound(%ld)=%g\n", i, upper_bound(i));
    }
	
	    
    return 0;
}


long maxlong(void){
	return LONG_MAX;
}


double upper_bound(long n){
	double fraction, power, UB;
	
	if(n>=6){
		fraction = n /2.0;
		power = pow(fraction, n);
		UB = power;
	}	 	 
	else
		UB = 719.0; 	
	
	return UB;
}
	
	
	
	
	

