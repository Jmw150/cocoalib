#ifndef CoCoA_PREPROCESSOR_DEFNS_H
#define CoCoA_PREPROCESSOR_DEFNS_H

// Created automatically by CoCoALib configure script.
// Created on date  2022-03-23  at time  15:39:35

// Some defines are used in src/AlgebraicCore/BuildInfo.C

// To enable debugging, reconfigure with option "--debug"
#undef CoCoA_DEBUG

// define CoCoA_MEMPOOL_DEBUG if you have memory problems
#undef CoCoA_MEMPOOL_DEBUG
// define CoCoA_MEMPOOL_DISABLE to avoid using MemPools
#undef CoCoA_MEMPOOL_DISABLE


#define CoCoA_ULONG2LONG 1
#define CoCoA_ULONGLONG2LONGLONG 1

#undef CoCoA_THREADSAFE_HACK
#define CoCoA_WITH_BOOST
#undef CoCoA_WITH_NORMALIZ
#undef CoCoA_WITH_FROBBY
#undef CoCoA_WITH_GFAN
#undef CoCoA_WITH_MATHSAT
#undef CoCoA_WITH_GSL

#endif
