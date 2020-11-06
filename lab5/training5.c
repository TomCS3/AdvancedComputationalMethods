#include <stdio.h>

/* function prototype */
long string_length(char s[]);

int main(void) {
  char s1[]="Hello";
  char s2[]="x";
  char s3[]="line 1\tline 2\n";

  printf("%20s | %s\n", "string_length(s)", "s");
  printf("%20ld | %s\n", string_length(s1), s1);
  printf("%20ld | %s\n", string_length(s2), s2);
  printf("%20ld | %s\n", string_length(s3), s3);
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
		

