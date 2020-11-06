#include <stdio.h>
#include <math.h>

#define N 10
#define XMAX 10
#define XMIN 1

int main(void){
	double x, f_x, i;
	i= (XMAX-XMIN)/(N-1.0);	   
	
	for(x=XMIN;x<=10;x+=i){
		f_x=sin(x);
		printf("%.6f %.6f\n", x, f_x);
	}

	return 0;
}

