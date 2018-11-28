/*

csvToFrags.c

Exports a CSV to frags

Use: csvToFrags file.csv file.frags

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#define  MAX_LEVELS  900000
#define CSV_BUFF 1024

int main(int ac,char** av){

	
	if(ac<3){
		printf("Use: csvToFrags <file.csv> <file.frags>\n");
		exit(-1);
	}
	
	// Read fragments
	struct FragFile frag;
	uint64_t xtotal,ytotal;
	//f=readFragmentsv2(av[1],&nf,&xtotal,&ytotal);


	FILE * fCSV = fopen(av[1], "rt");
	FILE * fFrags = fopen(av[2], "wb");
	if(fFrags == NULL){ fprintf(stderr, "Could not open output file\n"); exit(-1); }
	if(fCSV == NULL){ fprintf(stderr, "Could not open input file\n"); exit(-1); }


	char line[CSV_BUFF]; line[0] = '\0';
	while(strncmp(line, "Type", 4) != 0){
		fgets(line, CSV_BUFF, fCSV);
		if(strncmp(line, "SeqX length", 11) == 0){
			sscanf(line, "SeqX length          : %"PRIu64"\n", &xtotal);
		}  
		if(strncmp(line, "SeqY length", 11) == 0){
			sscanf(line, "SeqY length          : %"PRIu64"\n", &ytotal);
		}  
	}

	writeSequenceLength(&xtotal, fFrags);
	writeSequenceLength(&ytotal, fFrags);

	
	/*
	
	// Print header. Frags file info
	printf("All by-Identity Ungapped Fragments (Hits based approach)\n");
	printf("[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>\n");
	printf("SeqX filename        : Unknown\n");
	printf("SeqY filename        : Unknown\n");
	printf("SeqX name            : Unknown\n");
	printf("SeqY name            : Unknown\n");
	printf("SeqX length          : %"PRIu64"\n", xtotal);
	printf("SeqY length          : %"PRIu64"\n", ytotal);
	printf("Min.fragment.length  : 0\n");
	printf("Min.Identity         : 0\n");
	printf("Tot Hits (seeds)     : 0\n");
	printf("Tot Hits (seeds) used: 0\n");
	printf("Total fragments      : 0\n");
	printf("========================================================\n");
	printf("Total CSB: 0\n");
	printf("========================================================\n");
	printf("Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%%ident,SeqX,SeqY\n");

	*/
	float d1,d2;
				
	while(!feof(fCSV)){
		if(feof(fCSV)) break;
		
		fgets(line, CSV_BUFF, fCSV);
		sscanf(line, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%"PRIu64",%"PRIu64"\n",&frag.xStart, &frag.yStart, &frag.xEnd, &frag.yEnd, &frag.strand, &frag.block, &frag.length, &frag.score, &frag.ident, &d1, &d2, &frag.seqX, &frag.seqY);
		/*


		if(frag.strand=='r'){
			frag.yStart = ytotal - frag.yStart - 1;
			frag.yEnd = ytotal - frag.yEnd - 1;
		}

		*/
		
		
		writeFragment(&frag, fFrags);
				
	}

	fclose(fCSV);
	fclose(fFrags);
	
	return 0;
}
