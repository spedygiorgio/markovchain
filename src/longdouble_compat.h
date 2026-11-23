#ifndef MARKOVCHAIN_LONGDOUBLE_COMPAT_H
#define MARKOVCHAIN_LONGDOUBLE_COMPAT_H

/*
 * Portable accumulator type that follows R's long-double support.
 *
 * When R is built with --disable-long-double (CRAN's "noLD" check),
 * HAVE_LONG_DOUBLE is *not* defined in Rconfig.h.  On such builds R's
 * own sum() / colSums() accumulate in plain double, so package C++ code
 * must do the same to produce bit-identical results.
 *
 * Usage: replace every `long double` accumulator with `mc_ldbl`.
 */

#include <Rconfig.h>

#ifdef HAVE_LONG_DOUBLE
typedef long double mc_ldbl;
#else
typedef double mc_ldbl;
#endif

#endif /* MARKOVCHAIN_LONGDOUBLE_COMPAT_H */
