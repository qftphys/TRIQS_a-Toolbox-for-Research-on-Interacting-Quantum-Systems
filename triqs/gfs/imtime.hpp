/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2013 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./gf_classes.hpp"
#include "./meshes/matsubara_time.hpp"

namespace triqs {
namespace gfs {

 // singularity
 template <> struct gf_default_singularity<imtime, matrix_valued> {
  using type = tail;
 };
 template <> struct gf_default_singularity<imtime, scalar_valued> {
  using type = tail;
 };
 template <> struct gf_default_singularity<imtime, matrix_real_valued> {
  using type = tail;
 };
 template <> struct gf_default_singularity<imtime, scalar_real_valued> {
  using type = tail;
 };

 /// ---------------------------  HDF5 ---------------------------------
 template <typename Singularity> struct gf_h5_name<imtime, matrix_valued, Singularity> {
  static std::string invoke() { return "ImTime"; }
 };
 template <typename S, int R> struct gf_h5_name<imtime, tensor_valued<R>, S> : gf_h5_name<imtime, matrix_valued, S> {
  static std::string invoke() { return "ImTimeTv"+std::to_string(R); }
 };

 // If the data is real, write a real array, otherwise a complex array
 template <typename T, typename S, typename E> struct gf_h5_write_data<imtime, T, S, E> : gf_h5_write_data_real_or_complex_runtime{};

 /// ---------------------------  closest mesh point on the grid ---------------------------------

 template <typename Singularity, typename Target> struct gf_closest_point<imtime, Target, Singularity, void> {
  // index_t is int
  template <typename G, typename T> static int invoke(G const *g, closest_pt_wrap<T> const &p) {
   double x = double(p.value) + 0.5 * g->mesh().delta();
   int n = std::floor(x / g->mesh().delta());
   return n;
  }
 };

}
}

