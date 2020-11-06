#include <stdio.h>

#define MAXLINE 1000

long string_length(char s[]);
void lstrip(char s[]);


int main(void) {
  char test1[] = "   Hello World";

  printf("Original string reads  : |%s|\n", test1);
  lstrip(test1);
  printf("l-stripped string reads: |%s|\n", test1);

  return 0;
}

long string_length(char s[]){
	int i=0, len;
	while(s[i] != '\0'){
		i++;
	}
	len = i;
	return len;
}

void lstrip(char s[]){
	int len=string_length(s);
	int i=0, j=0;
	if (len>0 && len<=MAXLINE){
		while(s[i]=='\b' || s[i]=='\t' || s[i]==' '){
			i++;
		}
		for(;i<=len;j++, i++){
			s[j] = s[i];
		}
	}
} 

