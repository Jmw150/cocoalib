//   Copyright (c)  2017,2022  John Abbott and Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-resultant.H"

#include "CoCoA/ring.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/verbose.H"

//#include <vector>
using std::vector;
#include <iostream>
using std::endl;

namespace CoCoA
{

  namespace // anonymous
  {

    // Leading coeff viewed as poly in x
    RingElem lcx(ConstRefRingElem f, long x)
    {
      CoCoA_ASSERT(!IsZero(f));
      const vector<RingElem> c = CoeffVecWRT(f,indet(owner(f),x));
      return c.back();
    }

    // pseudo-remainder 
    RingElem prem(RingElem f, const RingElem& g, long x)
    {
      CoCoA_ASSERT(owner(f) == owner(g));
      CoCoA_ASSERT(!IsZero(g));
      if (IsZero(f)) return f;
      const ring& P = owner(f);
      const long degg = deg(g,x);
      if (degg == 0) return zero(P);
      long degf = deg(f,x);
      if (degf < degg) return f;
      const RingElem lcg = lcx(g,x);
      f *= power(lcg, degf-degg+1);
      while (degf >= degg)
      {
        const RingElem lcf = lcx(f,x);
        f -= (lcf/lcg)*g*power(indet(P,x), degf-degg); // exact division!
        if (IsZero(f)) return f;
        degf = deg(f,x);
      }
      return f;
    }

  } // end of namespace anonymous

  vector<RingElem> SubresultantSeq(ConstRefRingElem f, ConstRefRingElem g, long x)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NotPolyRing, "SubresultantSeq: args 1 & 2 must be in a poly ring");
    if (owner(g) != P)
      CoCoA_THROW_ERROR(ERR::MixedRings, "SubresultantSeq: args 1 & 2 must be in same ring");
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR(ERR::BadIndex, "SubresultantSeq: indet index out of range");
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR(ERR::NotNonZero, "SubresultantSeq: args 1 & 2 must be non-zero");
    VerboseLog VERBOSE("SubresultantSeq");
    VERBOSE(30) << "Inputs:" << endl;
    VERBOSE(30) << "f = " << f << endl;
    VERBOSE(30) << "g = " << g << endl;
    vector<RingElem> S;
    RingElem s = power(lcx(g,x), deg(f,x)-deg(g,x));
    RingElem A = g;
    RingElem B = prem(f,-g,x);
    if (IsZero(B)) { S.push_back(B); return S; } ////////// BUG BUG BUG:  DODGY HACK!!??!!
    VERBOSE(30) << "Before loop prem = " << B << endl;
    while (true)
    {
      if (IsZero(B)) return S;
      long dA = deg(A,x);
      long dB = deg(B,x);
      S.push_back(B);
      long delta = dA-dB;
      RingElem C;
      if (delta > 1)
      {
        VERBOSE(30) << "case delta > 1; in fact delta = " << delta << endl;
        VERBOSE(30) << "lcx(B,x) = " << lcx(B,x) << endl;
        VERBOSE(30) << "s = " << s << endl;
        VERBOSE(30) << "gcd = " << gcd(s, lcx(B,x)) << endl;
        C = power(lcx(B,x), delta-1)*B/power(s, delta-1);
        VERBOSE(30) << "C = " << C << endl;
        S.push_back(C);
      }
      else
      {
        VERBOSE(30) << "case delta=1" << endl;
        C = B;
      }
      if (dB == 0) return S;
      B = prem(A, -B, x);
      const RingElem D = power(s,delta)*lcx(A,x);
      VERBOSE(30) << "prem = " << B << endl;
      VERBOSE(30) << "denom = " << D << endl;
      B /= D; // exact div!
      A = C;
      s = lcx(A,x);
    }
    // never get here
  }


  RingElem resultant(ConstRefRingElem f, ConstRefRingElem g) // f,g univariate; result is in CoeffRing!
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NotPolyRing, "resultant: args 1 & 2 must be in a poly ring");
    if (owner(g) != P)
      CoCoA_THROW_ERROR(ERR::MixedRings, "resultant: args 1 & 2 must be in same ring");
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR(ERR::NotNonZero, "resultant: args 1 & 2 must be non-zero");
    const long x = UnivariateIndetIndex(f);
    if (x < 0 || (!IsConstant(g) && x != UnivariateIndetIndex(g)))
      CoCoA_THROW_ERROR(ERR::BadArg, "resultant: args 1 & 2 must be univariate polys");
    return ConstantCoeff(resultant(f,g,x));
  }

  RingElem resultant(ConstRefRingElem f, ConstRefRingElem g, long x)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NotPolyRing, "resultant: args 1 & 2 must be in a poly ring");
    if (owner(g) != P)
      CoCoA_THROW_ERROR(ERR::MixedRings, "resultant: args 1 & 2 must be in same ring");
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR(ERR::NotNonZero, "resultant: args 1 & 2 must be non-zero");
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR(ERR::BadIndex, "resultant: indet index out of range");

    return SubresultantSeq(f,g,x).back();
  }


  RingElem discriminant(ConstRefRingElem f)    ///< discr of univariate polynomial; result is in CoeffRing!
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NotPolyRing, "discriminant: arg 1 must be in a poly ring");
    if (IsConstant(f))
      CoCoA_THROW_ERROR(ERR::NotNonZero, "discriminant: arg 1 must be non-constant");
    const long x = UnivariateIndetIndex(f);
    if (x < 0)
      CoCoA_THROW_ERROR(ERR::BadArg, "discriminant: arg must be univariate");
    return ConstantCoeff(discriminant(f,x));
  }
  
  RingElem discriminant(ConstRefRingElem f, long x) ///< discr of multivariate polynomial
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::NotPolyRing, "resultant: arg 1 must be in a poly ring");
    if (IsZero(f))
      CoCoA_THROW_ERROR(ERR::NotNonZero, "discriminant: arg 1 must be non-zero");
    if (x < 0 || x >= NumIndets(P))
      CoCoA_THROW_ERROR(ERR::BadIndex, "discriminant: indet index out of range");

    const RingElem R = resultant(f,deriv(f,x),x);
    const long degf4 = deg(f)%4;
    const int sign = (degf4/2 == 0)?1:-1;
    return sign*(R/lcx(f,x));
  }

} // end of namespace CoCoA
  
