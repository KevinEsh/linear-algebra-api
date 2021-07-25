#ifndef LINEAR_ALGEBRA2_H
#define LINEAR_ALGEBRA2_H

#include <string>
#include <vector> 
#include <map>
#include <limits>
#include <iostream>
#include <fstream>
#include <cmath>
#include <regex>

class LinearAlgebra
{
private:
    int N1, N2, N3, solve_method, num_eigens;
    double* Matrix;
    double* inp_vector;
    double* down_diagonal;
    double* diagonal;
    double* up_diagonal;
    std::map<std::string, int> cache; //para evitar repetir calculos inesesarios
    bool pivot(int);
    void change_rows(int, int);
public:
    void change_ptr(double**,double**);
    double* out_vector;
    double* Q; 
    double* R;
    std::vector<double > eigenvalues;
    std::vector<double* > eigenvectors;
    LinearAlgebra(int, int);
    ~LinearAlgebra();
    void set_Matrix(std::string, bool);
    void elim_gausiana();
    void set_inp_vector(std::string, bool);
    void set_inp_vector(double*, int);
    void set_out_vector(double*, int);
    void set_Matrix_as_tridiagonal(double*, double*, double*);
    void solve();
    void solve_diag();
    void solve_triang_inf(int);
    void solve_triang_sup(int);
    void solve_tridiagonal_system();
    void print_array();
    void print_out_vector();
    void print_inp_vector();
    double L1_norm(double*, int);
    double Linf_norm(double*, int);
    double cuadratic_error(double*, int, double*, int);
    double relative_error(double*, int, double*, int);
    void sum_rows(double*);
    void sum_colms(double*);
    void LU_factorization();
    void LDL_factorization();
    void jacobbi(int, double);
    void gauss_seidel(int, double);
    bool is_diag_dom();
    void pow_method(double*, int, double);
    void inverse_pow_method(double*, int, double);
    void getall_eigens(std::string, int);
    void divide_for(double*, double, int);
    double dot(double*, double*, int);
    double max(double*, int);
    double max_abs(double*, int);
    double min(double*, int);
    double min_abs(double*, int);
    void norm_vec(double*, int);
    double* ones(int);
    void minus_proy(double*, double*, int); 
    void jacobi_eigens(double**, double**, bool, int);
    void rotation_jacobi_eigens(double*, double*, int, int, int, int);
    double max_abs_Mat(const double*, int, int, int&, int&, bool, bool);
    void QR_factorization();
    void subspace_eigens(int, int, double);
    double* dot_mat(const double*,int, int, const double*, int, int);
    double* transpose(double*, int, int);
    double frobenius_norm(double*, int, int, double);
    void gram_schmit(double*, int, int);
    void conjugated_gradient(double, double);
    double* copy(double*, int, int); 
    void add_mat(double*, int, int, double*, int, int, double*);
    void subtraction_mat(double*, int, int, double*, int, int, double*);
    void set_solve_method(std::string);
    void solve_for_inv_method();
};

#include "LinearAlgebra_def.h"

#endif