#include <stdio.h>
#include <stdlib.h>

void use_fib_array(long n);
long* make_long_array(long n);
long* make_fib_array(long n);

int main(void) {
	use_fib_array(10);
	return 0;
}

void use_fib_array(long N) {
  /* N is the maximum number for fibarray length */
  	long n;      /* counter for fibarray length */
  	long i;      /* counter for printing all elements of fibarray */
  	long *fibarray;  /* pointer to long -- pointer to the fibarray itself*/

	/* Print one line for each fibarray length n*/
  	for (n=2; n<=N; n++) {
    /* Obtain an array of longs with data */
    fibarray = make_fib_array(n);

    /* Print all elements in array */
    printf("fib(%2ld) : [",n);
    for (i=0; i<n; i++) {
      printf(" %ld", fibarray[i]);
    }
    printf(" ]\n");

    /* free array memory */
    free(fibarray);
  }
}

long* make_long_array(long n) {
	long *a;
	a = (long *)malloc(sizeof(long)*n);
	if (a == NULL) {
		printf("Memory allocation failed\n");
		return NULL;
	}
	else {
		printf("Allocated %ld bytes for p.\n", (long) sizeof(long)*n);
		return a;
		free(a);
	}
}
	
long* make_fib_array(long n){
	int i;
	if(make_long_array(n) == NULL){
		return NULL;
	}
	else{
		long *a = make_long_array(n);
		a[0] = 0;
		a[1] = 1;
		for(i=2;i<=n;i++){
			a[i] = a[i-1] + a[i-2];
		}  
		return a;
	}
}	 

