#include <stdio.h>

#define MAXLINE 1000

/* Function void rstrip(char s[])
modifies the string s: if at the end of the string s there are one or more spaces,
then remove these from the string.

The name rstrip stands for Right STRIP, trying to indicate that spaces at the 'right'
end of the string should be removed.
*/
long string_length(char s[]);
void rstrip(char s[]);


int main(void) {
  char test1[] = "";

  printf("Original string reads  : |%s|\n", test1);
  rstrip(test1);
  printf("r-stripped string reads: |%s|\n", test1);

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

void rstrip(char s[]){
	int len=string_length(s);
	int i=len-1;
	if (len>1 && len<=MAXLINE){
		while(s[i]=='\t' || s[i]=='\n' || s[i]==' '){
			i--;
		}
		s[++i] = '\0';
	}
} 

