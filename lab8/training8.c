#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void use_mix(void);
char* mix(char *s1, char *s2);


int main(void){
	use_mix();
	
	return 0;
}


void use_mix(void) {
    char s1[] = "Hello World";
    char s2[] = "1234567890!";

    printf("s1 = %s\n", s1);
    printf("s2 = %s\n", s2);
    printf("r  = %s\n", mix(s1, s2));
}

char* mix(char *s1, char *s2) {
	int len = strlen(s1);
	char *r;
	r = (char *)malloc(sizeof(char)*2*len);
	if(r == NULL) {
		return NULL;
	}
	else {
		int i, j;
		for(i=0, j=0; j<len; j++, i+=2) {
			r[i] = s1[j];
			r[i+1] = s2[j];
			printf("Value of i = %d\n", i);
			printf("Value of j = %d\n", j);
		}
		r[2*len] = '\0';
		printf("Value of len = %d\n", len);
		return r;
		
	}
}
		 
	

