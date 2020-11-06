#include <stdio.h>
#define MAXLINE 1000 /* maximum length of string */

/* function prototype */
long string_length(char s[]);
void reverse(char source[], char target[]);


int main(void) {
  char original[] = "This is a test: can you print me in reverse character order?";
  char reversed[MAXLINE];

  printf("%s\n", original);
  reverse(original, reversed);
  printf("%s\n", reversed);
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

/* reverse the order of characters in 'source', write to 'target'.
   Assume 'target' is big enough. */
void reverse(char source[], char target[]) {
	int i, j=0, len;
	len = string_length(source);
	if(len>0){
		for(i=len-1;i>=0;i--, j++){
			target[j] = source[i];
		} 
	}
	target[j]='\0';
}	   

