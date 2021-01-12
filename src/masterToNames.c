/*

masterToNames.c

Adds names/indices to the sequence IDs of a master csv file generated with CHROMEISTER/GECKO

Use: masterToNames master.csv indicesX indicesY

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
#define ID_LEN 50

uint64_t iterative_binary_search(uint64_t * array, int64_t start_index, int64_t end_index, uint64_t element){

	int64_t middle = 0;
	while (start_index <= end_index){
		middle = start_index + (end_index- start_index )/2;
		if (array[middle] == element) return (uint64_t) middle;
		if (array[middle] < element)  start_index = middle + 1;
		else end_index = middle - 1;
	}
	if(middle > end_index) return (uint64_t) middle - 1;
	if(element > array[middle]) return (uint64_t) middle;
	if(element < array[middle] && middle > 0) return (uint64_t) middle - 1;
	return (uint64_t) middle;
}

uint64_t asciiToUint64(const char *text)
{
    uint64_t number=0;

    for(int i=0; i<strlen(text); i++)
    {
        char digit=text[i]-'0';           
        number=(number*10)+(uint64_t)digit;
    }

    return number;
}

int main(int ac,char** av){

	
	if(ac<7){
		printf("Use: csvToFrags <file.csv> <indicesX> <indicesY> <seqs in X> <seqs in Y> <0 for num ID and 1 for names>\n");
		exit(-1);
	}
	
	// Read fragments
	struct FragFile frag;
	uint64_t xtotal,ytotal;
	//f=readFragmentsv2(av[1],&nf,&xtotal,&ytotal);


	FILE * fCSV = fopen(av[1], "rt");
	if(fCSV == NULL){ fprintf(stderr, "Could not open input file\n"); exit(-1); }

	FILE * fIX = fopen(av[2], "rt");
	if(fIX == NULL){ fprintf(stderr, "Could not open input X index file\n"); exit(-1); }
	FILE * fIY = fopen(av[3], "rt");
	if(fIY == NULL){ fprintf(stderr, "Could not open input Y index file\n"); exit(-1); }


	int seqsInX = atoi(av[4]);
	int seqsInY = atoi(av[5]);

	int useNames = atoi(av[6]);

	uint64_t indicesX[seqsInX];
	uint64_t indicesY[seqsInY];

	char seqNamesX[seqsInX][ID_LEN];
	char seqNamesY[seqsInY][ID_LEN];

	


	int pos = 0;
	char line[CSV_BUFF]; line[0] = '\0';
	char nameID[CSV_BUFF]; nameID[0] = '\0';
	char stringo[CSV_BUFF];
	while(!feof(fIX) && pos < seqsInX){
		fgets(line, CSV_BUFF, fIX);
		sscanf(line, "%s %s", &stringo, &nameID);
		//printf("%s %s %d\n", stringo, nameID, strlen(nameID));
		indicesX[pos] = asciiToUint64(stringo);
		strncpy(seqNamesX[pos], nameID, ID_LEN);
		seqNamesX[pos][ID_LEN-1] = '\0';
		++pos;
		nameID[0] = '\0';
	}
	fclose(fIX);

	line[0] = '\0'; pos = 0;
	nameID[0] = '\0';
    while(!feof(fIY) && pos < seqsInY){
        fgets(line, CSV_BUFF, fIY);
		sscanf(line, "%s %s", &stringo, &nameID);
        indicesY[pos] = asciiToUint64(stringo);
		strncpy(seqNamesY[pos], nameID, ID_LEN);
		seqNamesY[pos][ID_LEN-1] = '\0';
		++pos;
		nameID[0] = '\0';
    }
	fclose(fIY);


	line[0] = '\0';
	while(strncmp(line, "Type", 4) != 0){
		fgets(line, CSV_BUFF, fCSV);
		fprintf(stdout, "%s", line);

		if(strncmp(line, "SeqX length", 11) == 0){
			sscanf(line, "SeqX length          : %"PRIu64"\n", &xtotal);
		}  
		if(strncmp(line, "SeqY length", 11) == 0){
			sscanf(line, "SeqY length          : %"PRIu64"\n", &ytotal);
		}  
	}

	//writeSequenceLength(&xtotal, fFrags);
	//writeSequenceLength(&ytotal, fFrags);

	
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

		// Frag,75365314,1873657,75365476,1873819,f,0,163,484,142,74.23,0.87,0,0		
		fgets(line, CSV_BUFF, fCSV);
		sscanf(line, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%f,%f,%"PRIu64",%"PRIu64"\n",&frag.xStart, &frag.yStart, &frag.xEnd, &frag.yEnd, &frag.strand, &frag.block, &frag.length, &frag.score, &frag.ident, &d1, &d2, &frag.seqX, &frag.seqY);


		frag.seqX = iterative_binary_search(indicesX, 0, (int64_t) seqsInX - 1, frag.xStart);
		frag.seqY = iterative_binary_search(indicesY, 0, (int64_t) seqsInY - 1, frag.yStart);


		if(useNames == 0){
			fprintf(stdout, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%"PRIu64",%"PRIu64"\n",frag.xStart, frag.yStart, frag.xEnd, frag.yEnd, frag.strand, frag.block, frag.length, frag.score, frag.ident, d1, d2, frag.seqX, frag.seqY);
		}else{
			fprintf(stdout, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,%"PRId64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,%s,%s\n",frag.xStart, frag.yStart, frag.xEnd, frag.yEnd, frag.strand, frag.block, frag.length, frag.score, frag.ident, d1, d2, seqNamesX[frag.seqX], seqNamesY[frag.seqY]);
		}


		//printf(" for %"PRIu64" i found this: seqX:%"PRIu64"@ %"PRIu64"\n", frag.xStart, frag.seqX, indicesX[frag.seqX] );

		/*

		if(frag.strand=='r'){
			frag.yStart = ytotal - frag.yStart - 1;
			frag.yEnd = ytotal - frag.yEnd - 1;
		}

		*/
		
		
				
	}

	fclose(fCSV);
	
	return 0;
}
