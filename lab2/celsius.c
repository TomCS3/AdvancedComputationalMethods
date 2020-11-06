#include <stdio.h>

int main(void){
	int celsius;
	float fahrenheit;
	
	for(celsius = -30; celsius <= 30; celsius += 2){
		fahrenheit = celsius * 9.0 / 5.0 +32;
		printf("%3d = %5.1f\n", celsius, fahrenheit);
	}
	return 0;
}

