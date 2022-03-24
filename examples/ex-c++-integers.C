// Copyright (c) 2022  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;

#include <limits>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing the C++ integer type \"long\", which\n"
  "we recommend just for indexing; otherwise use \"BigInt\".      \n";

const string LongDescription =
  "This is an example showing the C++ integer type \"long\".  There  \n"
  "are other types with different bit-widths; all can be \"signed\"  \n"
  "or \"unsigned\" (just for non-negative values).  Computations     \n"
  "are fast, but values may overflow silently (on most computers). \n"
  "We give an example showing how awkward it can be to avoid or    \n"
  "anticipate overflow -- there are no such problem with BigInt!   \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // C++ offers many different integral types.  When working
  // with CoCoALib we recommend using just "long".
  // The type "char" is good for characters in strings.
  // We recommend avoiding values of "unsigned" integral types;
  // unfortunately several C++ functions return such values.
  // C++ does "automatic promotion" of integral values
  // [this is both good and bad].


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // We recommend using just the integral type "long";
    // this is actually an abbreviation for "signed long int".

    // Integer constants such as 3 or -99 are known as "literals"
    // (you may see this name in a compiler message).
    // For literals outside the range -9999 to 9999 it is a good
    // idea to append the suffix "L" -- a lower case "l" is also
    // allowed but looks confusingly similar to the digit "1".

    long NEGATIVE = -2147483647L; // OK on all binary computers
    long POSITIVE =  2147483647L; // OK on all binary computers
    cout << "All values between " << NEGATIVE << " and " << POSITIVE << " are OK" << endl;

    // Integers outside the range above may not work on some
    // computers.  The actual permitted range can be found like this:
    long MOST_NEGATIVE = numeric_limits<long>::min();
    long MOST_POSITIVE = numeric_limits<long>::max();
    if (MOST_NEGATIVE < NEGATIVE  ||  MOST_POSITIVE > POSITIVE)
      cout << "Actual range of representable values is wider." << endl << endl;

    // It is likely that -MOST_POSITIVE works as expected,
    // but that -MOST_NEGATIVE does not work as expected!


    // What are longs useful for?
    // (*) size of a container (e.g. vector or string)
    // (*) index into a container
    // (*) loops over elements in a container
    // (*) representing values from a limited range (e.g. years)

    // Arithmetic with longs:
    // (*) integer division (e.g.  2/3  gives  0)
    // (*) be careful that OVERFLOW does not occur
    //     [*YOU* must check -- the computer does not]
    // (*) if in doubt use BigInt (slower, but safer)
    
    // Silly example:  compute 1+2+3+...+n
    long n = 1000000L;  // Try 5000000000L as well...
    long sum = 0;
    for (long i=1; i <= n; ++i)
      sum += i;
    cout << "Sum of integers 1 up to " << n << "  is " << sum << endl;
    // sum will overflow if the result is greater than MOST_POSITIVE


    // --------------------------------------------
    // This next part is a bit trickier.
    // We want to compute the sum directly by using the formula:
    //   sum(1..n)  =  (1/2)*n*(n-1)
    long bad_formula    = n*(n+1)/2;
    long better_formula = (n*(n+1))/2;
    long best_formula   = (n%2 == 0) ? ((n+1)*(n/2)) : (n*((n+1)/2));

    cout << "bad_formula gives "    << bad_formula    << endl;
    cout << "better_formula gives " << better_formula << endl;
    cout << "best_formula gives "   << best_formula   << endl;

    // Before reading on: why did I call them as I did?

    // bad_formula:  it is unspecified when the division occurs:
    //     the compiler may choose to multiply first then divide,
    //     or it may choose to compute (n+1)/2 first then multiply
    //     [recall that (n+1)/2 implies integer division]
    
    // better_formula:  using brackets we instruct the compiler
    //      to muliply first, then divide; the problem is that
    //      the multiplication may overflow even if the answer
    //      is small enough to fit into a long.

    // best_formula:  "best" but not "most readable"!
    //      depending on whether n is even or odd we tell
    //      the compiler to divide first (but knowing that
    //      the division will be exact), and then multiply.
    //      This will overflow iff sum in the loop overflows.

    // !!CAUTION!!  Maybe bad_formula works...
    //   ...today
    //   ...on your current computer
    //   ...with your current compiler.
    // But it may give different results with a different computer/compiler.
    // It is not unambiguous portable code!  But it looks fine... argh!

    
    // --------------------------------------------
    // Take a deep breath...
    // [AIM: convince you to use BigInt for general integer computations]
    
    // How could we check if overflow will occur (or has occurred)?

    // Usually better to check beforehand, but the checks make the
    // code slower and harder to read.  Below is the code to evaluate
    // the above formula while checking for overflow (and printing out
    // a warning if overflow will occur).

    // For simplicity we ASSUME that n >= 0.
    long safest = 0;  // will be set properly by "if" below
    if (n == 0)
      safest = n; // handle the case n=0
    else if (n == MOST_POSITIVE)
        cerr << "!!! Overflow would occur !!!" << endl;
    // Henceforth we know that n+1 cannot cause overflow
    else if (n%2 == 0)
    {
      // Case n is even & n > 0
      long halfn = n/2;  // division is exact, cannot overflow, halfn > 0.
      if (MOST_POSITIVE/halfn >= n+1)  // division OK because halfn non-zero.
        safest = halfn*(n+1); // no overflow in (n+1); no overflow in product.
      else
        cerr << "!!! Overflow would occur !!!" << endl;
    }
    else
    {
      // Case n is odd, n > 0
      long halfn1 = (n+1)/2; // (n+1) cannot overflow, division is exact, halfn1 >= 0.
      if (MOST_POSITIVE/n >= halfn1) // division OK because n non-zero (but halfn1 may be zero)
        safest = halfn1*n; // no overflow in product.
      else
        cerr << "!!! Overflow would occur !!!" << endl;
    }
    // Correct answer is in "safest", or else we printed an overflow message.
    cout << "If no overflow message then correct answer is " << safest << endl;
    // Phew! (no, I did not get the code right on my first attempt)


    // ----- CONCLUSION -----
    //   It is best to use longs just for indexing;
    //   use BigInts for general arithmetic with integers.


    // For masochists:
    //  - look up integral types (char, short, int, long, long long)
    //  - look up integral promotions (on cppreference)
    //  - look up about unsigned integral types
  }

} // end of namespace CoCoA

// IGNORE THE STUFF BELOW (at least for now)

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-integers.C,v 1.3 2022/03/21 10:07:38 bigatti Exp $
// $Log: ex-c++-integers.C,v $
// Revision 1.3  2022/03/21 10:07:38  bigatti
// Summary: suggested  n = 5000000000L; exlicitely
//
// Revision 1.2  2022/03/17 18:29:18  abbott
// Summary: Better comments
//
// Revision 1.1  2022/03/17 16:16:01  abbott
// Summary: New examples
//
// Revision 1.7  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2021/10/07 14:27:45  abbott
// Summary: Added comment
//
// Revision 1.5  2021/09/03 12:58:41  abbott
// Summary: Improved comments
//
// Revision 1.4  2017/02/24 16:58:00  abbott
// Summary: Removed an empty line
//
// Revision 1.3  2017/02/15 12:21:53  abbott
// Summary: Added bool as a basic type
//
// Revision 1.2  2017/02/10 16:40:35  abbott
// Summary: Minor improvement
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
