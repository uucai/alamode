#pragma once

#include "pointers.h"
#include "kpoint.h"
#include <complex>
#include <Eigen/Dense>

namespace PHON_NS
{
    class DistList
    {
    public:
        unsigned int cell_s;
        double dist;

        DistList();

        DistList(const unsigned int cell_s_, const double dist_)
        {
            cell_s = cell_s_;
            dist = dist_;
        }
    };

    inline bool operator<(const DistList a, const DistList b)
    {
        return a.dist < b.dist;
    }

    struct ShiftCell
    {
    public:
        int sx, sy, sz;
    };

    struct MinimumDistList
    {
    public:
        double dist;
        std::vector<ShiftCell> shift;
    };

    struct KpointSymmetry
    {
    public:
        unsigned int symmetry_op;
        unsigned int knum_irred_orig;
        unsigned int knum_orig;
    };

    class Scph : protected Pointers
    {
    public:

        Scph(class PHON *phon);
        ~Scph();

        unsigned int kmesh_scph[3];
        unsigned int kmesh_interpolate[3];
        unsigned int ialgo;

        bool restart_scph;
        bool warmstart_scph;
        bool lower_temp;
        double tolerance_scph;

        double **xk_scph, **kvec_na_scph;
        double **xk_interpolate;
        std::vector<std::vector<KpointList>> kp_irred_scph;
        std::vector<std::vector<KpointList>> kp_irred_interpolate;

        void exec_scph();
        void setup_scph();
        void finish_scph();

        double mixalpha;
        unsigned int maxiter;
        bool print_self_consistent_fc2;
        bool selfenergy_offdiagonal;
        bool relax_coordinate;

    private:

        unsigned int nk_scph;
        unsigned int nk_interpolate;

        double ***vec_for_v3, *invmass_for_v3;
        double ***vec_for_v4, *invmass_for_v4;
        int **evec_index3;
        int **evec_index4;
        int *kmap_interpolate_to_scph;

        double **eval_harmonic;
        std::complex<double> ***evec_harmonic;
        std::complex<double> ****delta_dymat_scph;

        std::complex<double> im;
        std::complex<double> ***V4_inner;
        int ngroup, ngroup2;
        std::vector<double> *fcs_group;
        std::vector<double> *fcs_group2;
        unsigned int *knum_minus_scph;
        double ***eval_anharmonic;
        double **omega2_harmonic;
        std::complex<double> ****mat_transform_sym;
        std::vector<int> *small_group_at_k;
        std::vector<int> *symop_minus_at_k;
        KpointSymmetry *kpoint_map_symmetry;

        std::complex<double> *exp_phase, ***exp_phase3;
        int nk_grid[3];
        int nk_represent;
        unsigned int tune_type;
        double dnk[3];
        MinimumDistList ***mindist_list_scph;

        void setup_kmesh();
        void setup_eigvecs();
        void setup_pp_interaction();
        void setup_transform_ifc();
        void setup_transform_symmetry();
        void write_scph_energy(double ***);
        void write_scph_bands(double ***);
        void write_scph_dos(double ***);
        void write_scph_thermodynamics(double ***);
        void write_scph_msd(double ***, std::complex<double> ****);

        void load_scph_dymat_from_file(std::complex<double> ****);
        void store_scph_dymat_to_file(std::complex<double> ****);

        void exec_scph_main(std::complex<double> ****);
        void compute_V4_array_all(std::complex<double> ***,
                                  std::complex<double> ***, const bool, const bool);

        void compute_V4_array_all2(std::complex<double> ***,
                                   std::complex<double> ***, const bool);

        void compute_V3_array_all(std::complex<double> ***,
                                  std::complex<double> ***, const bool);

        void calc_new_dymat_with_evec(std::complex<double> ***,
                                      double **, std::complex<double> ***);

        void compute_anharmonic_frequency(std::complex<double> ***, double **,
                                          std::complex<double> ***, const double,
                                          std::vector<int> *, bool &, std::complex<double> ***,
                                          const bool);

        void exec_interpolation(std::complex<double> ***,
                                double **, std::complex<double> ***);
        void exec_interpolation2(std::complex<double> ***,
                                 double **, std::complex<double> ***);

        void r2q(double *, const unsigned int, const unsigned int,
                 const unsigned int, const unsigned int, const unsigned int,
                 std::complex<double> ***, std::complex<double> **);

        void diagonalize_interpolated_matrix(std::complex<double> **, double *,
                                             std::complex<double> **, const bool);

        void find_degeneracy(std::vector<int> *, const unsigned int,
                             std::vector<std::vector<KpointList>>, double **);

        double distance(double *, double *);
        void symmetrize_dynamical_matrix(const unsigned int, Eigen::MatrixXcd &);
        void replicate_dymat_for_all_kpoints(std::complex<double> ***);
        void duplicate_xk_boundary(double *, std::vector<std::vector<double>> &);

        void write_anharmonic_correction_fc2(std::complex<double> ****, const unsigned int);
        void mpi_bcast_complex(std::complex<double> ****,
                               const int, const int, const int);
    };

    extern "C"
    {
        void zgemm_(const char *transa, const char *transb, int *m, int *n,
                    int *k, std::complex<double> *alpha, std::complex<double> *a, int *lda,
                    std::complex<double> *b, int *ldb, std::complex<double> *beta,
                    std::complex<double> *c, int *ldc);

        void zgeev_(const char *jobvl, const char *jobvr, int *n, std::complex<double> *a, int *lda,
                    std::complex<double> *w, std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr,
                    std::complex<double> *work, int *lwork, double *rwork, int *info);
    }
}