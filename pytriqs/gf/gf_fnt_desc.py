from cpp2py.wrap_generator import *
import re

# This modules contains functions that may be called directly by users
m = module_(full_name = "pytriqs.gf.gf_fnt", doc = "C++ wrapping of functions on Green functions ...", app_name="triqs")

import meshes
import singularities

m.add_include("<triqs/gfs.hpp>")
m.add_include("<triqs/gfs/transform/pade.hpp>")
m.add_include("<triqs/gfs/legacy_for_python_api.hpp>")

m.add_include("<cpp2py/converters/vector.hpp>")
m.add_include("<triqs/cpp2py_converters.hpp>")

m.add_using("namespace triqs::arrays")
m.add_using("namespace triqs::gfs")
m.add_using("triqs::utility::mini_vector")
m.add_preamble("""
""")

# density
m.add_function("matrix<dcomplex> density(gf_view<imfreq, matrix_valued> g)",   doc = "Density, as a matrix, computed from a Matsubara sum")
m.add_function("matrix<dcomplex> density(gf_view<legendre, matrix_valued> g)", doc = "Density, as a matrix, computed from a Matsubara sum")
m.add_function("dcomplex density(gf_view<imfreq, scalar_valued> g)",   doc = "Density, as a complex, computed from a Matsubara sum")
m.add_function("dcomplex density(gf_view<legendre, scalar_valued> g)", doc = "Density, as a complex, computed from a Matsubara sum")


for Target in  ["scalar_valued", "matrix_valued", "tensor_valued<3>", "tensor_valued<4>"]:

    for Meshes in [["imtime", "imfreq"], ["retime", "refreq"], ["cyclic_lattice", "brillouin_zone"]]:
        # === Direct Fourier

        # Setter
        m.add_function("void set_from_fourier(gf_view<%s, %s> g_out, gf_view<%s, %s> g_in)"%(Meshes[1], Target, Meshes[0], Target),
                calling_pattern = "g_out = fourier(g_in)",
                doc = """Fills self with the Fourier transform of g_in""")
        # Factory function
        m.add_function(name = "make_gf_from_fourier",
                signature="gf_view<%s, %s> make_gf_from_fourier(gf_view<%s, %s> g_in)"%(Meshes[1], Target, Meshes[0], Target),
                doc ="""Create Green function from the Fourier transform of g_in""")

        # === Inverse Fourier

        # Setter
        m.add_function("void set_from_inverse_fourier(gf_view<%s, %s> g_out, gf_view<%s, %s> g_in)"%(Meshes[0], Target, Meshes[1], Target),
                calling_pattern = "g_out = inverse_fourier(g_in)",
                doc = """Fills self with the inverse Fourier transform of g_in""")
        # Factory function
        m.add_function(name = "make_gf_from_inverse_fourier",
                signature="gf_view<%s, %s> make_gf_from_inverse_fourier(gf_view<%s, %s> g_in)"%(Meshes[0], Target, Meshes[1], Target),
                doc ="""Create Green function from the inverse Fourier transform of g_in""")

    # make_real_in_tau
    m.add_function("gf_view<imfreq, %s> make_real_in_tau(gf_view<imfreq, %s> g)"%(Target, Target),
                doc = "Ensures that the Fourier transform of the Gf, in tau, is real, hence G(-i \omega_n)* =G(i \omega_n)")

    # is_gf_real_in_tau
    m.add_function("bool is_gf_real_in_tau(gf_view<imfreq, %s> g, double tolerance = 1.e-13)"%Target)


# fit_tail
m.add_function("void fit_tail(gf_view<imfreq, matrix_valued> g, tail_view known_moments, int max_moment, int n_min, int n_max, bool replace_by_fit = true)",
            doc = """Set the tail by fitting""")
m.add_function("void fit_tail(gf_view<imfreq, matrix_valued> g, tail_view known_moments, int max_moment, int neg_n_min, int neg_n_max, int pos_n_min, int pos_n_max, bool replace_by_fit = true)",
            doc = """Set the tail by fitting""")
m.add_function("void fit_tail(gf_view<imfreq, scalar_valued> g, __tail_view<scalar_valued> known_moments, int max_moment, int n_min, int n_max, bool replace_by_fit = true)",
            doc = """Set the tail by fitting""")
m.add_function("void fit_tail(gf_view<imfreq, scalar_valued> g, __tail_view<scalar_valued> known_moments, int max_moment, int neg_n_min, int neg_n_max, int pos_n_min, int pos_n_max, bool replace_by_fit = true)",
            doc = """Set the tail by fitting""")

# set_from_pade
m.add_function("void set_from_pade (gf_view<refreq, matrix_valued> gw, gf_view<imfreq, matrix_valued> giw, int n_points = 100, double freq_offset = 0.0)",
        calling_pattern = "pade(gw, giw, n_points, freq_offset)",
        doc = """""")

# set_from_legendre
m.add_function("void set_from_legendre(gf_view<imfreq, matrix_valued> gw, gf_view<legendre, matrix_valued> gl)",
            calling_pattern = "gw = legendre_to_imfreq(gl)",
            doc = """Fills self with the legendre transform of gl""")

m.add_function("void set_from_legendre(gf_view<imtime, matrix_valued> gt, gf_view<legendre, matrix_valued> gl)",
            calling_pattern = "gt = legendre_to_imtime(gl)",
            doc = """Fills self with the legendre transform of gl""")

# set_from_imfreq
m.add_function("void set_from_imfreq(gf_view<legendre, matrix_valued> gl, gf_view<imfreq, matrix_valued> gw)",
            calling_pattern = "gl = imfreq_to_legendre(gw)",
            doc = """Fills self with the legendre transform of gw""")

# set_from_imtime
m.add_function("void set_from_imtime(gf_view<legendre, matrix_valued> gl, gf_view<imtime, matrix_valued> gt)",
            calling_pattern = "gl = imtime_to_legendre(gt)",
            doc = """Fills self with the legendre transform of gt""")

# set_from_imfreq
m.add_function("void set_from_imfreq(gf_view<legendre, matrix_valued> gl, gf_view<imfreq, matrix_valued> gw)",
            calling_pattern = "gl = imfreq_to_legendre(gw)",
            doc = """Fills self with the legendre transform of gw""")

# enforce_discontinuity
m.add_function("void enforce_discontinuity(gf_view<legendre, matrix_valued> gl, matrix_view<double> disc)", doc = """Modify the coefficient to adjust discontinuity""")

# rebinning_tau
m.add_function("gf<imtime, matrix_valued> rebinning_tau(gf_view<imtime,matrix_valued> g, int new_n_tau)", doc = "Rebins the data of a GfImTime on a sparser mesh")

# GfLegendre specific functions
m.add_function("void enforce_discontinuity(gf_view<legendre, matrix_valued> gl, matrix_view<double> disc)", doc = """Modify the coefficient to adjust discontinuity""")

########################
##   Code generation
########################

m.generate_code()

