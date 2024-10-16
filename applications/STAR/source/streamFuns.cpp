#include "streamFuns.h"
#include "ErrorWarning.h"
#include <fstream>
#include <sys/statvfs.h>
#include <sys/stat.h>
#include <stdio.h>
#define fstream_Chunk_Max 2147483647
#include <cstring>

void createDirectory(const string dirPathIn, const mode_t dirPerm, const string dirParameter, Parameters &P)
{
	string dirPath = dirPathIn.substr(0,dirPathIn.find_last_of('/')+1);
	if (mkdir(dirPath.c_str(), dirPerm) == -1) {
		if ( errno == EEXIST ) {//directory exists
			P.inOut->logMain << dirParameter << " directory exists and will be overwritten: " << dirPath <<std::endl;
		} else {//other error
			//will try to create parent directories
			size_t i1=dirPath.find_first_of('/',1);
			while (i1<dirPath.size()) {
				string dirPath1=dirPath.substr(0,i1);
				if (mkdir(dirPath1.c_str(), dirPerm) == -1) {
					if ( !(errno == EEXIST) ) {//error
				        exitWithError("EXITING because of fatal OUTPUT FILE error: could not create output directory: " + dirPath1 +
				        			  " for " + dirParameter + " " + dirPathIn +
									  "\n ERROR: " + strerror(errno) +
				                      "\nSOLUTION: check the path and permissions.\n",
				                       std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
					};
				};
				i1=dirPath.find_first_of('/',i1+1);
			};
			P.inOut->logMain << dirParameter << " directory and its parents created: " << dirPath <<std::endl;
		};
    } else {
    	P.inOut->logMain << dirParameter << " directory created: " << dirPath <<std::endl;
    };
};

unsigned long long fstreamReadBig(std::ifstream &S, char* A, unsigned long long N) {
    unsigned long long C=0;
    for (unsigned long long ii=0; ii<N/fstream_Chunk_Max; ii++) {
        S.read(A+C,fstream_Chunk_Max);
        C+=S.gcount();
        if (!S.good()) break;
    };
    S.read(A+C,N%fstream_Chunk_Max);
    C+=S.gcount();
    return C;
};

void fstreamWriteBig(std::ofstream &S, char* A, unsigned long long N, std::string fileName, std::string errorID, Parameters &P) {

    struct statvfs statvfsBuf;
    statvfs(fileName.c_str(), &statvfsBuf);
    P.inOut->logMain << "Writing " << N << " bytes into " <<fileName << " ; empty space on disk = " << statvfsBuf.f_bavail * statvfsBuf.f_bsize <<" bytes ..." <<flush;

    unsigned long long C=0;
    unsigned long long iC;
    for (iC=0; iC<N/fstream_Chunk_Max; iC++) {
        S.write(A+C,fstream_Chunk_Max);
        C+=fstream_Chunk_Max;
    };
    if (!S.fail()) S.write(A+C,N%fstream_Chunk_Max);
    if (S.fail()) {//failed to write

        struct statvfs statvfsBuf;
        statvfs(fileName.c_str(), &statvfsBuf);

//         system(( "ls -lL "+ P.pGe.gDir + " > "+ P.pGe.gDir +"/error.info 2>&1").c_str());
//         ifstream error_info((P.pGe.gDir +"/error.info").c_str());
//         P.inOut->logMain <<error_info.rdbuf();

        struct stat statBuf;
        stat(fileName.c_str(), &statBuf);

        remove(fileName.c_str());

        ostringstream errOut;
        errOut << errorID<<": exiting because of *OUTPUT FILE* error: could not write the output file "<< fileName <<"\n";
        errOut << "fail()=" <<S.fail() <<" ; bad()="<< S.bad()<<"\n";
        errOut << "Error while trying to write chunk # " << iC << "; "<< C << " bytes\n";
        errOut << "File size full = "<< N <<" bytes\n";
        errOut << "File size on disk = " << statBuf.st_size<<" bytes\n";
        errOut << "Solution: check that you have enough space on the disk\n";
        errOut << "Empty space on disk = " << statvfsBuf.f_bavail * statvfsBuf.f_bsize <<" bytes\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_FILE_WRITE, P);
    };
    P.inOut->logMain << " done\n" <<flush;
};

std::ofstream &ofstrOpen (std::string fileName, std::string errorID, Parameters &P) {//open file 'fileName', generate error if cannot open
    std::ofstream & ofStream = *new std::ofstream(fileName.c_str(), std::fstream::out | std::fstream::trunc);
    if (ofStream.fail()) {//
        ostringstream errOut;
        errOut << errorID<<": exiting because of *OUTPUT FILE* error: could not create output file "<< fileName <<"\n";
        errOut << "SOLUTION: check that the path exists and you have write permission for this file. Also check ""ulimit -n"" and increase it to allow more open files.\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
    };
    return ofStream;
};

std::fstream &fstrOpen (std::string fileName, std::string errorID, Parameters &P, bool flagDelete) {//open file 'fileName', generate error if cannot open
    //std::fstream &fStream = *new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out );
    //std::fstream &fStream = *new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);

    std::fstream *fStreamP;
    if (flagDelete) {//truncate the file if it exists
        fStreamP = new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);
    } else {//try to open exising file
        fStreamP=new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out );
        if (fStreamP->fail()) //did not work <= file does not exist => open with trunc (the above command does not work on new file)
            fStreamP = new std::fstream(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);
    };
    
    if (fStreamP->fail()) {//
        ostringstream errOut;
        errOut << errorID<<": exiting because of *OUTPUT FILE* error: could not create input/output file "<< fileName <<"\n";
        errOut << "Solution: check that the path exists and you have write permission for this file\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
    };
    return *fStreamP;
};

std::ifstream & ifstrOpen (std::string fileName, std::string errorID, std::string solutionString, Parameters &P) {
    //open file 'fileName', generate error if cannot open
    std::ifstream & ifStream = *new std::ifstream(fileName.c_str());
    if (ifStream.fail()) {//
        ostringstream errOut;
        errOut << errorID<<": exiting because of *INPUT FILE* error: could not open input file "<< fileName <<"\n";
        errOut << "Solution: check that the file exists and you have read permission for this file\n";
        if (solutionString.size()>0) {
            errOut << "          "<< solutionString <<"\n";
        };
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
    };
    return ifStream;
};

ifstream & ifstrOpenGenomeFile (std::string fileName, std::string errorID, Parameters &P) {
     //open one of the genome files
     return ifstrOpen(P.pGe.gDir+"/"+fileName, errorID,  "if this file is missing from the genome directory, you will need to *re-generate the genome*", P);
};

void copyFile(string fileIn, string fileOut)
{//copy fileIn into FileOut
    std::ifstream  src(fileIn, std::ios::binary);
    std::ofstream  dst(fileOut,   std::ios::binary);
    dst << src.rdbuf();
};
