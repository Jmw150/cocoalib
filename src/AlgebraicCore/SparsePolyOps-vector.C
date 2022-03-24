//   Copyright (c)  2020  John Abbott,  Anna Bigatti
//   Original author: 2020 Julian Danner (interreduced transcoded from CoCoA-5)

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/SparsePolyOps-vector.H"

#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/geobucket.H" // for myMul
#include "CoCoA/interrupt.H"
#include "CoCoA/verbose.H"
#include "CoCoA/ReductionCog.H"


//#include <vector>
using std::vector;


namespace CoCoA
{

  // Deliberately not defined!
  // void interreduce(std::vector<RingElem>& v)
  // {
  //   std::swap(v, interreduced(v));
  // }


  // naive impl (orig transcoded from CoCoA-5 by Julian Danner)
  std::vector<RingElem> interreduced(std::vector<RingElem> v)
  {
    static const char* const FnName = "interreduced";
    if (v.empty()) { return v; } // ??? or error???
    if (!HasUniqueOwner(v)) CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    VerboseLog VERBOSE(FnName);
    //delete possible zeros in v
    const ring& P = owner(v[0]);
    v.erase(std::remove(v.begin(), v.end(), zero(P)), v.end());

    // this local fn is used in call to sort
    const auto CompareLPPs = [](const RingElem& f, const RingElem& g) { return LPP(f)<LPP(g); };
    long count = 0;
    while (true)
    {
      VERBOSE(90) << "round " << ++count << std::endl; // NB *always* incrs count!
      sort(v.begin(), v.end(), CompareLPPs);
      vector<RingElem> ans;
      RingElem rem;
      bool NewLPPfound = false;
      for (const auto& f: v)
      {
        CheckForInterrupt(FnName);
        rem = NR(f, ans);
        if (IsZero(rem)) continue;
        ans.push_back(rem);
        NewLPPfound = (NewLPPfound || LPP(rem) != LPP(f));
//        if (!NewLPPfound && LPP(rem) != LPP(f))
//          NewLPPfound = true;
      }
      if (!NewLPPfound) return ans;
      swap(v, ans); // quicker than: v = ans;
    }
  }


  RingElem IndetsProd(const std::vector<RingElem>& L)
  {
    if (L.empty()) CoCoA_THROW_ERROR(ERR::Empty, "IndetsProd");
    const ring& P = owner(L[0]);
    if (!IsSparsePolyRing(P)) CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "IndetsProd");
    const vector<long> VarIndices = IndetsIn(L);
    RingElem ans = one(P);
    for (auto k: VarIndices)
      ans *= indet(P,k);
    return ans;
  }


  // indices of indets appearing in (non-empty) L
  std::vector<long> IndetsIn(const std::vector<RingElem>& L)
  {
    if (L.empty()) CoCoA_THROW_ERROR(ERR::Empty, "IndetsIn");
    if (!IsSparsePolyRing(owner(L[0]))) CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "IndetsIn");
    const SparsePolyRing& P = owner(L[0]);
    const int nvars = NumIndets(P);
    const long n = len(L);
    for (long i=1; i < n; ++i)
      if (owner(L[i]) != P) CoCoA_THROW_ERROR(ERR::MixedRings, "IndetsIn");
    vector<bool> IndetSeen(nvars);
    int NumSeen = 0;
    for (long i=0; i < n; ++i)
    {
      const vector<long> IndetsInThisPoly = IndetsIn(L[i]);
      for (long xj: IndetsInThisPoly)
      {
        if (IndetSeen[xj]) continue;
        IndetSeen[xj] = true;
        if (++NumSeen == nvars) goto double_break;
      }
    }
    double_break:
    vector<long> ans; ans.reserve(nvars); // potentially wasteful reserve
    for (int i=0; i < nvars; ++i)
      if (IndetSeen[i]) ans.push_back(i);
    return ans;
  }


  //----------------------------------------------------------------------
  //??? the following functions to compute NR will be replaced by GBMill

  int FindReducerIndex(ConstRefPPMonoidElem pp, const std::vector<RingElem>& v)
  {
    const long nelems = len(v);
    for (long i=0; i < nelems; ++i)
      if (IsDivisible(pp, LPP(v[i])))
        return i;
    return -1;
  }


  inline int FindReducerIndex(const ReductionCog& F, const std::vector<RingElem>& v)
  {
    if ( IsActiveZero(F) ) return -1;
    return FindReducerIndex(ActiveLPP(F), v);
  }


  void ReduceActiveLM(ReductionCog& F, const std::vector<RingElem>& v)
  {
    int i;
    while ( (i = FindReducerIndex(F, v) ) != -1)
    {
      CheckForInterrupt("ReduceActiveLM");
      F->myReduce(v[i]);
    }
  }


  void reduce(ReductionCog& F, const std::vector<RingElem>& v)
  {
    ReduceActiveLM(F, v);
    while ( !IsActiveZero(F) )
    {
      F->myMoveToNextLM();
      ReduceActiveLM(F, v);
    }
  }
  //--------------------------------------------

  RingElem NR(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    const char* const FnName = "NR";
    const ring& P = owner(f);
    if (!IsPolyRing(P) || !IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR("Must be in polyring over a field", FnName);
    if (!HasUniqueOwner(v) || (!v.empty() && owner(v.front()) != P))
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if ( IsZero(f) ) return f;
    if ( v.empty() ) return f;
    RingElem ans(f);
    ReductionCog F = NewRedCogGeobucketField(owner(ans));
    F->myAssignReset(ans);
    reduce(F, v);
    F->myRelease(ans);
    return ans;
  }

  // naive implementation
  namespace { // anonymous

    RingElem DivAlgLM(vector<RingElem>& QuotRem, ConstRefRingElem f, const vector<RingElem>& v)
    {
      const SparsePolyRing P = owner(f);
      RingElem m(P);
      RingElem r =f;
      
      long j = FindReducerIndex(LPP(r), v);
      while (j != -1)
      {
        P->myDivLM(raw(m), raw(r), raw(v[j])); // m =LM(r)/LM(v[i]); no checks
        QuotRem[j] += m;
        r -= m * v[j];
        if (IsZero(r)) break;
        j = FindReducerIndex(LPP(r), v);
      }
      return r;
    }

  } // anonymous
  
  // RingElem NormalRemainder(ConstRefRingElem f, const vector<RingElem>& v)
  // {
  //   if (IsZero(f)) return f;
  //   const SparsePolyRing P = owner(f);
  //   RingElem ansNR(P);
  //   RingElem tmpNR =f;
  
  //   tmpNR = NRLM(f, v);
  //   while (!IsZero(tmpNR))
  //   {
  //     P->myMoveLMToBack(raw(ansNR), raw(tmpNR));
  //     tmpNR = NRLM(tmpNR, v);
  //   }
  //   return ansNR;
  // }


  std::vector<RingElem> TmpDivAlg(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    const char* const FnName = "NR";
    const SparsePolyRing& P = owner(f);
    if (!IsPolyRing(P) || !IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR("Must be in polyring over a field", FnName);
    if (!HasUniqueOwner(v) || (!v.empty() && owner(v.front()) != P))
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    //    if ( IsZero(f) ) return f;
    //    if ( v.empty() ) return f;
    RingElem ansNR(P);
    RingElem tmpNR =f;
    vector<RingElem> QuotRem(len(v)+1, zero(P));
    while (!IsZero(tmpNR))
    {
      tmpNR = DivAlgLM(QuotRem, tmpNR, v);
      if (IsZero(tmpNR)) break;
      P->myMoveLMToBack(raw(ansNR), raw(tmpNR)); // ansNR+=LM(tmpNR)
    }
    swap(QuotRem[len(v)], ansNR);
    return QuotRem;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-vector.C,v 1.3 2022/03/09 07:57:52 bigatti Exp $
// $Log: SparsePolyOps-vector.C,v $
// Revision 1.3  2022/03/09 07:57:52  bigatti
// Summary: added TmpDivAlg
//
// Revision 1.2  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2022/02/14 14:47:04  bigatti
// Summary: interreduced moved into SparsePolyOps-vector
//
// Revision 1.3  2020/11/19 18:29:10  abbott
// Summary: Added static keyword
//
// Revision 1.2  2020/10/23 07:54:20  abbott
// Summary: Improved commented out code (?!?)
//
// Revision 1.1  2020/10/14 20:01:54  abbott
// Summary: Renamed SparsePolyOps-interreduce to SparsePolyOps-interreduced
//
// Revision 1.2  2020/10/05 19:24:54  abbott
// Summary: Added comment
//
// Revision 1.1  2020/10/02 19:04:55  abbott
// Summary: New fn interreduce
//
//
//
