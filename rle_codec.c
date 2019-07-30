/*
	Donald Elmore
	RLE Codec
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compress(char* file) {
	unsigned char runCount = 0;
	unsigned char A, B;

	//open file for reading
	FILE *original = fopen(file, "rb");
	FILE *compressed = fopen("compressed", "w");

	if (original == NULL) {
		printf("Unable to read %s\n", file);
		exit(0);
	}

	//compression
	fread(&A, 1, 1, original);
	while(fread(&B, 1, 1, original) == 1) {
		if (A != B) {
			fwrite(&runCount, sizeof(runCount), 1, compressed);
			fwrite(&A, sizeof(A), 1, compressed);
			A = B;
			runCount = 0;
		} else if (runCount >= 255) {
			fwrite(&runCount, sizeof(runCount), 1, compressed);
			fwrite(&A, sizeof(A), 1, compressed);
			A = B;
			runCount = 0;
		} else {
			runCount++;
		}
	}

	fwrite(&runCount, sizeof(runCount), 1, compressed);
	fwrite(&A, sizeof(A), 1, compressed);

	fclose(original);
	fclose(compressed);
}

void decompress(char* file) {
	unsigned char A, B;
	int count, i;

	//open file for reading
	FILE *original = fopen(file, "rb");
	FILE *decomp = fopen("decompressed", "w");

	if (original == NULL) {
		printf("Unable to read %s\n", file);
		exit(0);
	}

	//decompression
	while(fread(&A, 1, 1, original) == 1) {
		fread(&B, 1, 1, original);
		count = (int)A;
		for (i = 0; i < count+1; i++) {
			fwrite(&B, sizeof(B), 1, decomp);
		}
	}
	fclose(original);
	fclose(decomp);
}

int main(int argc, char *argv[]) {
	//check that command line is correct
	if (argc != 3){
		printf("Error! Formatting is: \"executable\" [filename] [c/d]\n");
		printf("Where c is for compression, d is for decompression\n");
		exit(0);
	}
	if (argc != 3 && ((strcmp(argv[2],"c") == 0) || (strcmp(argv[2],"d") == 0))) {
		printf("Error! Formatting is: \"executable\" [filename] [c/d]\n");
		printf("Where c is for compression, d is for decompression\n");
		exit(0);
	}
	if (strcmp(argv[2],"c") == 0) 
		compress(argv[1]);
	if (strcmp(argv[2],"d") == 0)
		decompress(argv[1]);

	return 0;
}
