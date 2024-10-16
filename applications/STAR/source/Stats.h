#ifndef STATS_DEF
#define STATS_DEF

#include "IncludeDefine.h"
#include "Transcript.h"
#include "Parameters.h"


class Stats {
    public:
        uint readN;//number of reads from the file
        uint readBases;//number of input bases
//         uint mateLmax[2], mateLmin[2];//mates' max and min length

        uint mappedReadsU, mappedReadsM;
        uint mappedBases, mappedMismatchesN, mappedInsN, mappedDelN, mappedInsL, mappedDelL;
        double mappedPortion; //portion of the read length that has been mapped

        uint splicesN[SJ_MOTIF_SIZE];//non-can,3*can,annotated
        uint splicesNsjdb;

        uint unmappedOther, unmappedShort, unmappedMismatch, unmappedMulti, unmappedAll;

        uint chimericAll;

        time_t timeStart, timeStartMap, timeFinishMap, timeLastReport, timeFinish;
        
        Stats ();
        void resetN();
        void printShort(ostream*);
        void transcriptStats(Transcript &T, uint Lread);
        void addStats(Stats &S);
        void progressReportHeader(ofstream &progressStream);
        void progressReport(ofstream &progressStream) ;
        void reportFinal(ofstream &streamOut);
        void writeLines(ofstream &streamOut, const vector<int> outType, const string commStr, const string outStr);// write commented lines to text files with stats
        
        void qualHistCalc(const uint64 imate, const char* qual, const uint64 len);
        //void qualHistCalcSolo(const uint64 imate, const char* qual, const vector<uint32> stlen);
};
#endif
