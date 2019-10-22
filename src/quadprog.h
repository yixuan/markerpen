/* Copyright (C) 2010 */
/* Berwin A Turlach <Berwin.Turlach@gmail.com> */

/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 2 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along */
/* with this program; if not, write to the Free Software Foundation, Inc., */
/* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

#ifndef BAT_QUADPROG_H
#define BAT_QUADPROG_H

#include <R.h>

extern "C" {

void F77_SUB(qpgen2)
     (double *dmat, double *dvec, int *fddmat, int *n,
      double *sol, double *lagr, double *crval,
      double *amat, double *bvec, int *fdamat, int *q,
      int *meq, int *iact, int *nact, int *iter,
      double *work, int *ierr);
}

#endif
