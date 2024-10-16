#include "BAMbinSortByCoordinate.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "BAMfunctions.h"
#include "SequenceFuns.h"

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &genome, Solo &solo) {

    if (binS==0) return; //nothing to do for empty bins
    //allocate arrays
    char *bamIn=new char[binS+1];
    uint *startPos=new uint[binN*3];

    uint bamInBytes=0;
    //load all aligns
    for (uint it=0; it<nThreads; it++) {
        string bamInFile=dirBAMsort+to_string(it)+"/"+to_string((uint) iBin);
        ifstream bamInStream;
        bamInStream.open(bamInFile.c_str(),std::ios::binary | std::ios::ate);//open at the end to get file size
        int64 s1=bamInStream.tellg();
        if (s1>0)         {
            bamInStream.seekg(std::ios::beg);
            bamInStream.read(bamIn+bamInBytes,s1);//read the whole file
        } else if (s1<0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: failed reading from temporary file: " << dirBAMsort+to_string(it)+"/"+to_string((uint) iBin);
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        bamInBytes += bamInStream.gcount();
        bamInStream.close();
        remove(bamInFile.c_str());
    };
    if (bamInBytes!=binS) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
        errOut << "Expected bin size=" <<binS <<" ; size on disk="<< bamInBytes <<" ; bin number="<< iBin <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
    };

    //extract coordinates

    for (uint ib=0,ia=0;ia<binN;ia++) {
        uint32 *bamIn32=(uint32*) (bamIn+ib);
        startPos[ia*3]  =( ((uint) bamIn32[1]) << 32) | ( (uint)bamIn32[2] );
        startPos[ia*3+2]=ib;
        ib+=bamIn32[0]+sizeof(uint32);//note that size of the BAM record does not include the size record itself
        startPos[ia*3+1]=*( (uint*) (bamIn+ib) ); //read order
        ib+=sizeof(uint);
    };

    //sort
    qsort((void*) startPos, binN, sizeof(uint)*3, funCompareArrays<uint,3>);

    BGZF *bgzfBin;
    bgzfBin=bgzf_open((dirBAMsort+"/b"+to_string((uint) iBin)).c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
    if (bgzfBin==NULL) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal ERROR: could not open temporary bam file: " << dirBAMsort+"/b"+to_string((uint) iBin) << "\n";
        errOut <<"SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };

    outBAMwriteHeader(bgzfBin,P.samHeaderSortedCoord,genome.chrNameAll,genome.chrLengthAll);
    //send ordered aligns to bgzf one-by-one
    char bam1[BAM_ATTR_MaxSize];//temp array
    for (uint ia=0;ia<binN;ia++) {
        char* bam0=bamIn+startPos[ia*3+2];
        uint32 size0=*((uint32*) bam0)+sizeof(uint32);
        
        if (solo.pSolo.samAttrYes)
            solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->addBAMtags(bam0,size0,bam1);
        
        bgzf_write(bgzfBin, bam0, size0);
    };

    bgzf_flush(bgzfBin);
    bgzf_close(bgzfBin);
    //release memory
    delete [] bamIn;
    delete [] startPos;
};
