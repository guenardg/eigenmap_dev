/*************************************************************************
 
 (c) 2008-2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Registering routines and dynamic symbols**
 
 This file is part of eigenmap
 
 eigenmap is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 eigenmap is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with eigenmap.  If not, see <https://www.gnu.org/licenses/>.
 
 *************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void dist_geo_hvs(double*, double*, int*, int*, double*, double*);
extern void dist_geo_vif(double*, double*, int*, int*, double*, int*,
                         double*, double*, int*, double*);
extern void dist_Euclid(double*, double*, int*, int*, int*, double*, int*);

static const R_CMethodDef CEntries[] = {
  {"dist_geo_hvs",  (DL_FUNC) &dist_geo_hvs,  6},
  {"dist_geo_vif",  (DL_FUNC) &dist_geo_vif, 10},
  {"dist_Euclid",   (DL_FUNC) &dist_Euclid,   7},
  {NULL,            NULL,                     0}
};

static const R_CallMethodDef CallEntries[] = {
  {NULL,            NULL,                     0}
};

void R_init_eigenmap(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
