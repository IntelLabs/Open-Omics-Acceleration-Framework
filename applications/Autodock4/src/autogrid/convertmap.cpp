#include <stdio.h>
char shorten_num(double);

main(int argc, char *argv[])
{
    double numin;
    FILE *infile, *outfile;
    int ix;
    char buffer[80];

    infile = fopen(argv[1], "r");
    outfile = fopen(argv[2], "w");

    for(ix = 1; ix <= 6; ix++) {
	fgets(buffer, 81, infile);
	fputs(buffer, outfile);
    }

    while(fscanf(infile, "%lf", &numin) != EOF) {
	fprintf(outfile, "%c\n", shorten_num(numin));
    }

    fclose(infile);
    fclose(outfile);

    return 0;
}

char shorten_num(double num) 
{
    /****  MANIPULATION CODE  ****/
    char numout;

    if (num == 0.){
	numout = 0;
    } else if ((-12.8 < num) && (num < 0.)) {
	numout = num * 10.;
    } else if ((0. < num) && (num < 1280.)) {
	numout = num / 10.;
    } else if (num >= 1280.)
	numout = 127;
    else
	numout = -128;

    return numout;
}
