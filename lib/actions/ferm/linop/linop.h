// -*- C++ -*-
// $Id: linop.h,v 1.16 2004-01-12 14:54:08 bjoo Exp $

/*! \file
 * \brief Linear operators
 *
 * Various fermion linear operators
 */

/*! \defgroup linop Linear operators
 * \ingroup actions
 *
 * Various fermion linear operators
 */

#ifndef __linop_h__
#define __linop_h__

#include "lmdagm.h"
#include "lopscl.h"

#ifdef CHROMA_BUILD_WILSON
#include "linop_w.h"
#elif CHROMA_BUILD_STAGGERED
#include "linop_s.h"
#endif

#endif


