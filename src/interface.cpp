//    Copyright 2021 Jose Bouza
//
//    This file is part of Persistence Landscape Toolbox (PLT).
//
//    PLT is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published
//    by the Free Software Foundation, either version 2.1 of the License, or (at
//    your option) any later version.
//
//    PLT is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with PLT.  If not, see <http://www.gnu.org/licenses/>.
#include "pl.h"
#include <Rcpp.h>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(PersistenceLandscapeInterface)
RCPP_MODULE(Landscape) {

  class_<PersistenceLandscapeInterface>("PersistenceLandscape")
      .constructor<NumericMatrix, bool, double, double, double, double>()

      .method("getExact",
              &PersistenceLandscapeInterface::getPersistenceLandscapeExact,
              "Returns the PL in the exact representation.")
      .method("getDiscrete",
              &PersistenceLandscapeInterface::getPersistenceLandscapeDiscrete,
              "Returns the PL in the discrete representation.")
      .method("getInternal", &PersistenceLandscapeInterface::getInternal,
              "Returns the internal tensor representation of the PL.")
      .method("add", &PersistenceLandscapeInterface::sum,
              "Adds this PL to another.")
      .method("scale", &PersistenceLandscapeInterface::scale,
              "Scales this PL by a scaler.")
      .method("inner", &PersistenceLandscapeInterface::inner,
              "Take inner product of this PL with another.");

  Rcpp::function("PLsum", &PLsum);
  Rcpp::function("PLscale", &PLscale);
  Rcpp::function("PLinner", &PLinner);
  Rcpp::function("PLaverage", &PLaverage);
}

