#! /bin/bash

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# Auxiliary script for CoCoALib configure script.
# This script expects the env variables CXX and CXXFLAGS to be set.

# This script tries to check that CXX is a working C++ compiler,
# and that CXXFLAGS are suitable flags for it.  The check entails
# compiling a very simple source file (with and without CXXFLAGS).


if [ $# -ne 0 ]
then
  echo "ERROR: expected no args.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Check environment variable CXX
if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Check that CXX is an executable file:
FULLCXX=`which "$CXX" 2>/dev/null`
if [ $? -ne 0 -o \! \( -x "$FULLCXX" -a -r "$FULLCXX" -a -f "$FULLCXX" \) ]
then
  echo "ERROR: specified compiler \"$CXX\" is not executable.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir verify-compiler`

pushd "$TMP_DIR" > /dev/null

# Here is the simple source code we shall use to test the compiler:
/bin/cat > test-compiler-gnu.C <<EOF
#include <iostream>
using namespace std;
int main()
{
#ifdef __GNUC__
  cout << "gnu";
#else
  cout << "not gnu";
#endif
}
EOF

# Try plain compiler (without CXXFLAGS):
$CXX test-compiler-gnu.C -o test-compiler-gnu  > LogFile  2>&1
if [ $? -ne 0 -o \! -f test-compiler-gnu -o \! -x test-compiler-gnu ]
then
  echo "ERROR: Are you sure \"$CXX\" is a C++ compiler?   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi
/bin/rm test-compiler-gnu  # not necessary, just being tidy :-)

# Try compiler with CXXFLAGS:
$CXX $CXXFLAGS test-compiler-gnu.C -o test-compiler-gnu  > LogFile  2>&1
if [ $? -ne 0 -o \! -f test-compiler-gnu -o \! -x test-compiler-gnu ]
then
  echo "ERROR: Compilation flags \"$CXXFLAGS\" seem to be unsuitable for \"$CXX\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

COMPILER_TYPE=`./test-compiler-gnu`

# Clean up TMP_DIR
popd > /dev/null
/bin/rm -rf "$TMP_DIR"

echo "$COMPILER_TYPE"
