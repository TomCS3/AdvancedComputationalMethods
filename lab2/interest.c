#include <stdio.h>

int main(void){
	double s = 1000, debt = s, rate = 0.03, interest, total_interest=0, frac;
	int month;
	
	for(month = 1; month<=24; month++){
		interest = debt * rate;
		debt += interest;
		total_interest=debt-1000;
		frac=total_interest/10.0;
		printf("month %2d: debt=%7.2f, interest=%4.2f, total_interest=%7.2f, frac=%6.2f%%\n", month, debt, interest, total_interest, frac);
	}
	
	return 0;

}

