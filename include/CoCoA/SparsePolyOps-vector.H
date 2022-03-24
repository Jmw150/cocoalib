#ifndef CoCoA_SparsePolyOps_vector_H
#define CoCoA_SparsePolyOps_vector_H

//   Copyright (c)  2020-2021  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/ring.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  std::vector<RingElem> interreduced(std::vector<RingElem> v);
// NOT DEFINED!!  void interreduce(std::vector<RingElem>& v);

  std::vector<RingElem> TmpDivAlg(ConstRefRingElem f, const std::vector<RingElem>& v);
  // the remainder is in result[len(v)]   (this interface might change)

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-vector.H,v 1.3 2022/03/09 07:57:07 bigatti Exp $
// $Log: SparsePolyOps-vector.H,v $
// Revision 1.3  2022/03/09 07:57:07  bigatti
// Summary: added TmpDivAlg
//
// Revision 1.2  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2022/02/14 14:57:31  bigatti
// Summary: moved interreduced into SparsePolyOps-vector
//
// Revision 1.1  2020/10/14 20:01:54  abbott
// Summary: Renamed SparsePolyOps-interreduce to SparsePolyOps-interreduced
//
// Revision 1.1  2020/10/09 11:49:10  abbott
// Summary: New file for interreduced
//
//
//

#endif
