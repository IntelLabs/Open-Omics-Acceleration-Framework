#include "Variation.h"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "TimeFunctions.h"
#include "serviceFuns.cpp"
#include "ErrorWarning.h"

Variation::Variation (Parameters &Pin, vector <uint> &chrStartIn, map <string,uint> &chrNameIndexIn, bool yesVar) : P(Pin), chrStart(chrStartIn), chrNameIndex(chrNameIndexIn) {
    if (yesVar) {
        yes=true;

        //not used yet
        //varOutFileName=P.outFileNamePrefix+"Variation.out";
        //varOutStream.open(varOutFileName);

        vcfFile=P.var.vcfFile;
        loadVCF(vcfFile);
    } else{
        yes=false;
    };
};

void scanVCF(ifstream& vcf, Parameters& P, SNP& snp, vector <uint> &chrStart, map <string,uint> &chrNameIndex) {
    snp.N=0;
    uint nlines=0;
    uint32 nHomoz=0;//number of homoz alleles
    while (true) {
        string chr,id, ref, alt, dummy, sample;
        uint pos;
        nlines++;

        vcf >> chr;
        if (!vcf.good()) {
            break;
        };

        if (chr.at(0)!='#') {
            vcf >> pos >> id >> ref >> alt >> dummy >> dummy >> dummy >> dummy >> sample;

            vector <string> altV(3);

            if (ref.size()==1 && splitString(alt,',',altV)==1) {//only SNVs allowed - ref=1-char, alt could be comma separated list of 1-char. splitString returns the max lenght of the split strings
                altV.insert(altV.begin(),ref);//add ref to the beginning

                if (chrNameIndex.count(chr)==0) {//chr not in Genome
                    P.inOut->logMain << "WARNING: while processing varVCFfile file=" << P.var.vcfFile <<": chromosome '"<<chr<<"' not found in Genome fasta file\n";
                } else if (sample.size()<3) {
                    //undefined genotype
                } else if (sample.size()>3 && sample.at(3)!=':') {
                    P.inOut->logMain << "WARNING: while processing varVCFfile file=" << P.var.vcfFile <<": more than 2 alleles per sample for line number "<<nlines<<"\n";
                } else if (sample.at(0)=='0' && sample.at(2)=='0') {
                    //both alleles are reference, no need to do anything
                } else if (altV.at( atoi(&sample.at(0)) ).at(0)==ref.at(0) && altV.at( atoi(&sample.at(2)) ).at(0)==ref.at(0)) {
                    //both alleles are reference, no need to do anything
                    //this is a strange case in VCF when ALT allele(s) are equal to REF
                } else if ( P.var.heteroOnly && sample.at(0)==sample.at(2) )  {
                    nHomoz++;
                } else {
                    array<char,3> nt1;
                    nt1[0]=convertNt01234( ref.at(0) );
                    nt1[1]=convertNt01234( altV.at( atoi(&sample.at(0)) ).at(0) );
                    nt1[2]=convertNt01234( altV.at( atoi(&sample.at(2)) ).at(0) );
                    if (nt1[0]<4 && nt1[1]<4 && nt1[2]<4) {//only record if variant is ACGT
                        snp.lociV.push_back(pos-1+chrStart[chrNameIndex[chr]]);
                        snp.nt.push_back(nt1);
                        snp.N++;
                    };
                };
            };
        };
        getline(vcf,dummy);
    };

    if (nHomoz>0) {
        P.inOut->logMain << "WARNING: while processing varVCFfile file=" << P.var.vcfFile 
                         <<": found " << nHomoz << " homozygous SNPs, they will NOT be used." <<endl;
    };
};


void Variation::loadVCF(string fileIn) {
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... loading variations VCF\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ..... loading variations VCF\n" <<flush;

    ifstream & vcf = ifstrOpen(fileIn, ERROR_OUT, "SOLUTION: check the path and permissions of the VCF file: "+fileIn, P);
    scanVCF(vcf, P, snp, chrStart, chrNameIndex);

    snp.loci=new uint[snp.N];
    for (uint ii=0;ii<snp.N;ii++)
        snp.loci[ii]=snp.lociV[ii];
    snp.lociV.clear();


    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ..... Loaded VCF data, found "<<snp.N<< " SNPs"<<endl;
    if (snp.N==0) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL INPUT FILE ERROR: could not find any SNPs in VCF file: " <<fileIn<< "\n";
        errOut <<"SOLUTION: check formatting of the VCF file; unzip VCF file or use process substitution.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };


    uint *s1=new uint[2*snp.N];
    for (uint ii=0;ii<snp.N; ii++) {
        s1[2*ii]=snp.loci[ii];
        s1[2*ii+1]=ii;
    };

    qsort((void*)s1, snp.N, 2*sizeof(uint), funCompareUint1);

    vector< array<char,3> > nt1=snp.nt;
    for (uint ii=0;ii<snp.N; ii++) {
        snp.loci[ii]=s1[2*ii];
        snp.nt[ii]=nt1.at(s1[2*ii+1]);
    };
    //sort SNPs by coordinate
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ..... Finished sorting VCF data"<<endl;
};

void SNP::snpOnBlocks(uint blockStart, uint blockL, int blockShift, vector<vector<array<int,2>>> &snpV) {
    int32 isnp=binarySearch1b <uint> (blockStart, loci, N);
    while ((uint)isnp<N && loci[isnp]<(blockStart+blockL)) {
        for (int ii=0;ii<2;ii++) {
            if (nt[isnp][ii+1]!=nt[isnp][0]) {//allele different from reference
                array<int,2> snp1;
                snp1[0]=(int) (loci[isnp]-blockStart)+blockShift;
                snp1[1]=(int) nt[isnp][ii+1];
                snpV[ii].push_back(snp1);
            };
        };
        ++isnp;
    };
};

vector<vector<array<int,2>>> Variation::sjdbSnp(uint sjStart, uint sjEnd, uint sjdbOverhang1) {
    vector<vector<array<int,2>>> snpV(2);

    if (!yes) {//no variation, return 1 empty element
        vector<vector<array<int,2>>> snpV1(1);
        return snpV1;
    };

    snp.snpOnBlocks(sjStart-sjdbOverhang1, sjdbOverhang1, 0,            snpV);
    snp.snpOnBlocks(sjEnd+1,              sjdbOverhang1, sjdbOverhang1, snpV);

    if (snpV.at(0).empty() && snpV.at(1).empty()) {
        snpV.pop_back();
    } else if (snpV.at(0) == snpV.at(1)) {
        snpV.pop_back();
    };

    return snpV;
};
