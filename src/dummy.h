


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2020 -- 2020 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#ifndef rfutils_dummy_H
#define rfutils_dummy_H 1

typedef unsigned int Uint;

#ifdef DO_PARALLEL
#define HAS_PARALLEL true
#else
#define HAS_PARALLEL false
#endif
#if defined AVX2
#define HAS_AVX2 true
#else
#define HAS_AVX2 false
#endif
#if defined AVX
#define HAS_AVX true
#else
#define HAS_AVX false
#endif
#if defined SSSE3
#define HAS_SSSE3 true
#else
#define HAS_SSSE3 false
#endif
#if defined SSE2
#define HAS_SSE2 true
#else
#define HAS_SSE2 false
#endif
#if defined SSE
#define HAS_SSE true
#else
#define HAS_SSE false
#endif



#define AttachMessageN 1000
#define HAS_ONE_RELEVANT (HAS_PARALLEL || (HAS_AVX2 && AVX2_USED) || (HAS_AVX && AVX_USED) || (HAS_SSSE3 && SSSE3_USED) ||  HAS_SSE2 || HAS_SSE)
#define HAS_ALL_RELEVANT (HAS_PARALLEL && (HAS_AVX2 || !AVX2_USED) && (HAS_AVX || !AVX_USED) && (HAS_SSSE3 || !SSSE3_USED) &&  HAS_SSE2 && HAS_SSE)

#ifdef WIN32
#define AttachMessageX(PKG, HINT, AND)					\
  "'"#PKG"' %.20s %.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s.%.320s%.120s", \
    HAS_ONE_RELEVANT ? "sees" : "does not see any of",			\
    HAS_PARALLEL ? "OMP" : "",					\
    HAS_AVX2 && AVX2_USED ? ", AVX2" : "",				\
    HAS_AVX && AVX_USED ? ", AVX" : "",					\
    HAS_SSSE3 && SSSE3_USED ? ", SSSE3" : "",				\
    HAS_SSE2 && SSE2_USED ? ", SSE2" : "",				\
    HAS_SSE && SSE_USED ? ", SSE" : "",					\
    !HAS_ONE_RELEVANT || HAS_ALL_RELEVANT ? "" : ", but not ",		\
    !HAS_PARALLEL ? "OMP, " : "",					\
    !HAS_AVX2 && AVX2_USED ? "AVX2" : "",				\
    !HAS_AVX && AVX_USED ? ", AVX" : "",				\
    !HAS_SSSE3 && SSSE3_USED ? ", SSSE3" : "",				\
    !HAS_SSE2 && SSE2_USED ? ", SSE2" : "",				\
    !HAS_SSE && SSE_USED ? ", SSE" : "",				\
    HINT && ((!HAS_AVX2 && AVX2_USED) || (!HAS_AVX && AVX_USED) || !HAS_SSE2) ? "\nBy default '"#PKG"' is compiled with flag '-mavx' under your OS.\nIf you are sure that AVX2 is available, consider adding the flag '-march=native'\nto 'PKG_CXXFLAGS' in the file src/Makefile.win and then recompile\n'"#PKG"' "#AND"." : "", \
    HINT && (!HAS_AVX2 && AVX2_USED) ?					\
    "\nOr: try adding flag '-mavx2' to 'PKG_CXXFLAGS'" : ""
#else
#define AttachMessageX(PKG, HINT, AND) 				\
  "'"#PKG"' %.20s %.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s.%.320s%.120s%.120s%.350s", \
    HAS_ONE_RELEVANT ? "sees" : "does not see any of",			\
    HAS_PARALLEL ? "OMP" : "",						\
    HAS_AVX2 && AVX2_USED ? ", AVX2" : "",				\
    HAS_AVX && AVX_USED ? ", AVX" : "",					\
    HAS_SSSE3 && SSSE3_USED ? ", SSSE3" : "",				\
    HAS_SSE2 && SSE2_USED ? ", SSE2" : "",				\
    HAS_SSE && SSE_USED ? ", SSE" : "",					\
    !HAS_ONE_RELEVANT || HAS_ALL_RELEVANT ? "" : ", but not ",		\
    !HAS_PARALLEL ? "OMP, " : "",					\
    !HAS_AVX2 && AVX2_USED ? "AVX2" : "",				\
    !HAS_AVX && AVX_USED ? ", AVX" : "",				\
    !HAS_SSSE3 && SSSE3_USED ? ", SSSE3" : "",				\
    !HAS_SSE2 && SSE2_USED ? ", SSE2" : "",				\
    !HAS_SSE && SSE_USED ? ", SSE" : "",				\
    HINT && ((!HAS_AVX2 && AVX2_USED) || (!HAS_AVX && AVX_USED) || !HAS_SSE2) ? "\nWithout appropriate SIMD instruction set, the calculations might be slow.\nConsider recompiling '"#PKG"' "#AND" with flags e.g.,\n install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS=-march=native\")" : "", \
    HINT && (!HAS_AVX2 && AVX2_USED) ?					\
    "\n install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS=-mavx2\")" \
    : "",								\
    HINT && (!HAS_AVX && AVX_USED) ?					\
    "\n install.packages(\""#PKG"\" , configure.args=\"CXX_FLAGS=-mavx\")"\
    : "",								\
    HINT && ((!HAS_AVX2 && AVX2_USED) || (!HAS_AVX && AVX_USED) || !HAS_SSE2) ? "\nAlternatively consider installing '"#PKG"'\nfrom https://github.com/schlather/"#PKG", i.e.,\n   install.packages(\"devtools\")\n   library(devtools)\n   devtools::install_github(\"schlather/"#PKG"/pkg\")" : ""
#endif

#if defined ownprefixN
#define AttachMessage(PKG, HINT)  AttachMessageX(PKG, HINT, )
#else
#define AttachMessage(PKG, HINT)  AttachMessageX(PKG, HINT, )
#endif
  
#define ReturnAttachMessage(PKG,HINT) 	\
  SEXP Ans = PROTECT(allocVector(STRSXP, 1));	\
  char simd_msg[AttachMessageN];			\
  SPRINTF(simd_msg, AttachMessage(PKG,HINT)); \
  SET_STRING_ELT(Ans, 0, mkChar(simd_msg));		\
  UNPROTECT(1);						\
  return Ans;

#endif
