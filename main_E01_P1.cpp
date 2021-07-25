#include <iostream>
#include <fstream>
#include <string>
#include "LinearAlgebra.h"

using namespace std;

int main(){

    LinearAlgebra *M = new LinearAlgebra(1000, 1000);
    M->set_Matrix("A_1000.mtx");
    ofstream M_inv("Ainv_1000.mtx");

    //Metodo de la inversa con factorizacion LU
    M->LU_factorization();
    double* ei = new double[1000];
    for(int i=0; i <1000; i++) ei[i] = 0;

    for(int i =0; i<1000; i++){
        ei[i] = 1.0;
        M->set_inp_vector(ei, 1000);
        M->solve();
        
        for (int j=0; j < 1000; j++) M_inv<<to_string(M->out_vector[i])+" ";
        M_inv<<endl;

        ei[i] = 0.0;
    }
    delete M;

    LinearAlgebra *O = new LinearAlgebra(1000,1000);
    O->set_Matrix("Ainv_1000.mtx");
    O->set_inp_vector("b_1000.vec");
    O->set_solve_method("inv_method");
    O->solve();
    O->print_out_vector();
    return 0;
}