#ifndef CODE_BAMbinSortUnmapped
#define CODE_BAMbinSortUnmapped
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Genome.h"
#include "Solo.h"

#include SAMTOOLS_BGZF_H

void BAMbinSortUnmapped(uint32 iBin, uint nThreads, string dirBAMsort, Parameters &P, Genome &genome, Solo &solo);

#endif
