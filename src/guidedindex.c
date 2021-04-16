/********************* 


USAGE		<sequences.fasta>	The fasta sequences
			<index>		The output index name (text)

 * ------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "commonFunctions.h"
#include "structs.h"



int main(int ac, char **av) {

	FILE   *f, *index;

	if (ac!=3) terror("USE: guidedindex <sequences.fasta> <output>");
	
	if ((f=fopen64(av[1],"rt"))==NULL) terror("Could not open sequences fasta file");

	if ((index=fopen64(av[2],"wt"))==NULL) terror("Could not open output file");

	char c;
	uint64_t pos = 0, l = 0;
	while(!feof(f)){
        c = getc(f);
        if(c == '>'){
			if(l > 0) fprintf(index, " %"PRIu64"\n", l);
            fprintf(index, "%"PRIu64" ", pos);
            c = getc(f);
            while(c != '\n')
            {
				if(c == ' ') c = '_';
				if(c == ',') c = '_';
                fprintf(index, "%c", c);
                c = getc(f);
            }
			l = 0;
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            ++pos; ++l;
        }
    }

	fprintf(index, " %"PRIu64"\n", l);


	fclose(f);
	fclose(index);


	return 0;
}



