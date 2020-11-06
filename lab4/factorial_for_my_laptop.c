#include <stdio.h>
#include <limits.h>
#include <math.h>


/*function prototypes*/
long maxlong(void);
long factorial(long n);
double upper_bound(long n);

int main(void) {
    long i;

    /* The next line should compile once "maxlong" is defined. */
    printf("maxlong()=%ld\n", maxlong());

    /* The next code block should compile once "upper_bound" is defined. */

    
    for (i=-1; i<22; i++) {
        printf("upper_bound(%ld)=%g\n", i, upper_bound(i));
		printf("factorial(%ld)=%ld\n\n", i, factorial(i));
    }
	
	    
    return 0;
}


long maxlong(void){
	return LONG_MAX;
}


long factorial(long n){
	long j = (n-1); 
	double fac = n, ans;
	if(n<0)
		fac = -2;
	else if(n==0)
		fac = 1;
	else if(n==1)
		fac = 1;
	/*else if(n>19) passes with 100% because it tests as if memory on uni computers
		fac=-1;*/ 
	else{
		for(;j>0; j--){
			fac*=j;
		}
	}
	printf("n = %ld", n);

	if(fac > upper_bound(n)){
		fac = -1;
	}
	ans = (long)fac;
	return ans;
	
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
	if(UB > LONG_MAX)
		UB = LONG_MAX; 	   
	
	return UB;
}
	
	
	
	
	

