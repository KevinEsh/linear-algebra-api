using namespace std;

LinearAlgebra::LinearAlgebra(int n1, int n2)
{   
    N1 = n1, N2 = n2, N3 = n2; 
    Matrix = nullptr, inp_vector = nullptr, out_vector = nullptr;

    //Alocacion de memoria.
    try{
    Matrix = new double[N1*N2];
    inp_vector = new double[N3];
    out_vector = (double*)malloc(N3*sizeof(double));
    }
    catch(bad_alloc&){
        cout<<"Problemas en la alocación de memoria. Todas las funciones serán desactivadas"<<endl;
        solve_method =-1;
        return;
    }
    solve_method = 1; //se resolverá el sistema con elim_gaussiana por defecto. 
    cache["diag_dom"] = 0; // 0: aun no se sabe, 1: false, 2: true
}

LinearAlgebra::~LinearAlgebra()
{
    delete[] Matrix;
    delete[] inp_vector;
    delete[] out_vector;
    delete[] Q;
    delete[] R;
    for (int i=0; i<num_eigens; i++) delete[] eigenvectors[i];
}

void LinearAlgebra::set_Matrix(string file_name, bool header=false){
    ifstream inp_file (file_name);
    string line;
    double value;
    int i=0, j=0;
    regex pattern_of_numbers("[+-]*[0-9]+\\.?[0-9]*");

    if (!inp_file.is_open()){
        cout<<"Archivo "<<file_name<<" no encontrado. Intente de nuevo"<<endl;
        return;
    }
    if (header){
        getline(inp_file, line);
    }

    while(getline(inp_file, line)){
        regex_iterator<string::iterator> numbers ( line.begin(), line.end(), pattern_of_numbers);
        regex_iterator<string::iterator> endline;
        while (numbers!=endline){
            istringstream string_to_double(numbers->str());
            string_to_double >> value;
            Matrix[i*N2+j] = value; 
            j++;
            ++numbers;
        }
        i++;
        j=0;
    }
    inp_file.close();
}

void LinearAlgebra::set_inp_vector(string file_name, bool header=false){
    ifstream inp_file(file_name);
    string line;
    double value;
    int i=0;
    regex pattern_of_numbers("[+-]*[0-9]+\\.?[0-9]*");

    if (!inp_file.is_open()){
        cout<<"Archivo "<<file_name<<" no encontrado. Intente de nuevo"<<endl;
        return;
    }
    if (header){
        getline(inp_file, line);
    }
    while(getline(inp_file, line)){
        regex_iterator<string::iterator> numbers ( line.begin(), line.end(), pattern_of_numbers);
        regex_iterator<string::iterator> endline;
        istringstream string_to_double(numbers->str());
        string_to_double >> value;
        inp_vector[i] = value; 
        i++;
    }
    inp_file.close();    
}

void LinearAlgebra::set_inp_vector(double* vector, int size){
    if (N3 < size){
        free(out_vector);
        delete[] inp_vector;
        N3 = size;
        inp_vector = new double[N3];
    }
    for (int i=0; i<N3; i++) inp_vector[i] = vector[i]; //Lo pasamos a la memoria heap
}

void LinearAlgebra::set_out_vector(double* vector, int n3){
    if (N3 < n3){
        delete[] out_vector;
        N3 = n3;
        out_vector = new double[N3];
    }
    for (int i=0; i<N3; i++) out_vector[i] = vector[i]; //Lo pasamos a la memoria heap
}

void LinearAlgebra::solve_diag(){
    if (N1!=N2){
        printf("Matriz no es cuadrada, en cambio dimension (%d,%d) fue dada",N1,N2);
        return;
    }
    else if (N1!=N3){
        printf("Matriz (%d,%d) y vector (%d) no son compatibles",N1,N2,N3);
        return;
    }

    double deno;
    for (int i=0; i<N1; i++){
        deno = Matrix[i*N2+i];
        if (deno == 0.0){
            cout<<"El sistema de ecuaciones no tiene solucion, columna de zeros fue encontrada"<<endl;
            return;
        }
        out_vector[i] = inp_vector[i]/deno;
    }
}
void LinearAlgebra::solve_triang_sup(int diag = 0)
{   
    if (N1!=N2){
        printf("Matriz no es cuadrada, en cambio dimension (%d,%d) fue dada",N1,N2);
        return;
    }
    else if (N1!=N3){
        printf("Matriz (%d,%d) y vector (%d) no son compatibles",N1,N2,N3);
        return;
    }
    double sum, deno;
    //Reduccion hacia atras
    if (diag==0) deno = Matrix[(N1-1)*N2 + N1-1];
    else deno = 1;
    if (deno == 0.0){
        cout<<"El sistema de ecuaciones no tiene solucion, columna de zeros fue encontrada"<<endl;
        return;
    }
    out_vector[N1-1] = inp_vector[N1-1]/deno;
    for (int i = N1-2; i >=0 ; i--){
        if (diag==0) deno = Matrix[i*N2+i];
        else deno=1;
        if (deno == 0.0){
            cout<<"El sistema de ecuaciones no tiene solucion, columna de zeros fue encontrada"<<endl;
            return;
        }
        sum=0.0;   
        for (int j = i+1; j < N2; j++){
            sum += Matrix[i*N2 + j]*out_vector[j];        
        }
        out_vector[i] = (inp_vector[i] - sum)/deno;
    }
}

void LinearAlgebra::elim_gausiana()
{   
    if (N1!=N2){
        printf("Matriz no es cuadrada, en cambio dimension (%d,%d) fue dada",N1,N2);
        return;
    }
    if (N1!=N3) { 
        printf("Numero de filas de matriz y vector deben de ser iguales. En cambio (%d,%d) y (%d) fue dado", N1,N2,N3);
        return;
    }
    double deno;
    double nume;
    
    for (int k = 0; k < N1-1;k++){
        pivot(k);
        deno = Matrix[k*N2+k];
        if (deno == 0.0){
            if (!pivot(k)){//si no encontro un pivote diferente de cero, regresa
                cout<<"No se puede resolver el sistema de ecuaciones"<<endl;
                return;
            }
            deno = Matrix[k*N2+k]; //denominador actualizado
        }

        for (int i = k+1; i < N1; i++){
            nume = Matrix[i*N2+k];
            for (int j = k; j < N2; j++){ 
                Matrix[i*N2+j] -= (nume*Matrix[k*N2+j])/deno;
            }    
            inp_vector[i] -= (nume*inp_vector[k])/deno; 
        }
    }
    solve_triang_sup();
}

void LinearAlgebra::solve_triang_inf(int diag = 0)
{   
    if (N1!=N2){
        printf("Matriz no es cuadrada, en cambio una de dimension (%d,%d) fue dada",N1,N2);
        return;
    }
    else if (N1!=N3){
        printf("Matriz A (%d,%d) y vector b (%d) no son compatibles",N1,N2,N3);
        return;
    }
    
    double sum, deno;
    //Reduccion hacia atras
    if (diag==0) deno = Matrix[0];
    else deno = 1;
    if (deno == 0.0){
        cout<<"El sistema de ecuaciones no tiene solucion, columna de zeros fue encontrada"<<endl;
        return;
    }
    out_vector[0] = inp_vector[0]/deno;
    for (int i = 1; i < N1; i++){   
        if (diag==0) deno = Matrix[i*N2+i];
        else diag = 1;
        if (deno == 0.0){
            cout<<"El sistema de ecuaciones no tiene solucion, columna de zeros fue encontrada"<<endl;
            return;
        }
        sum = 0.0; 
        for (int j = 0; j < i; j++){
            sum += Matrix[i*N2+j]*out_vector[j];        
        }
        out_vector[i] = (inp_vector[i] - sum)/deno;
    }
}

void LinearAlgebra::print_array()
{   
    for (int i=0; i<N1; i++){
        cout<<"[";
        for (int j=0; j<N2; j++){
            cout<<Matrix[i*N2+j]<<", ";
        }
        cout<<"]"<<endl;
    }
    cout<<endl;
}

void LinearAlgebra::print_inp_vector()
{   
    cout<<"[";
    for (int i=0; i<N3; i++) cout<<inp_vector[i]<<", ";
    cout<<"]"<<endl<<endl;
}

void LinearAlgebra::print_out_vector()
{   
    cout<<"[";
    for (int i=0; i<N3; i++) cout<<out_vector[i]<<", ";
    cout<<"]"<<endl<<endl;
};

void LinearAlgebra::change_rows(int row1, int row2){
    double aux;
    for (int k = 0; k < N2; k++){
        aux = Matrix[row1*N2+k];
        Matrix[row1*N2+k] = Matrix[row2*N2+k];
        Matrix[row2*N2+k] = aux;

        aux = inp_vector[row1];
        inp_vector[row1] = inp_vector[row2];
        inp_vector[row2] = aux; 
    }
}

bool LinearAlgebra::pivot(int row1){
    double max = 0, elem;
    int row2, zeros=0;
    for (int i = row1; i<N1; i++){
        elem = abs(Matrix[i*N2+row1]); //barre la fila desde la diagonal
        if (elem == 0.0) zeros++; 
        else if (elem > max){
            max = elem;
            row2 = i;
        }
    }
    if (N1-row1==zeros){
        cout<<"Matriz singular, fila de zeros fue encontrada."<<endl;
        return false; //no hay fila sustituible
    }
    change_rows(row1,row2);
    return true;
}

double LinearAlgebra::cuadratic_error(double* vector1, int n1, double* vector2, int n2){
	if (n1!=n2){
		printf("Vectores no son de la misma dimension, (%d) y (%d) fueron dados\n",n1,n2);
		return 0.0;
	}
	double error=0;
	for (int i=0; i<n1; i++) error+= pow(vector1[i]-vector2[i],2);
	return error;    
}
double LinearAlgebra::relative_error(double* vector1, int n1, double* vector2, int n2){
	if (n1!=n2){
		printf("Vectores no son de la misma dimension, (%d) y (%d) fueron dados\n",n1,n2);
		return 0.0;
	}
	double error=0;
	for (int i=0; i<n1; i++) error+= abs(vector1[i]-vector2[i])/abs(vector1[i]);
	return error;
}	

double LinearAlgebra::L1_norm(double* inp_vector1, int n1){
	double norm=0;
	for (int i=0; i<n1; i++){
		norm+= abs(inp_vector[i]);
	}
	return norm;
}

double LinearAlgebra::Linf_norm(double* inp_vector1, int n1){
	double norm = 0, elem;
	for (int i=0; i<n1; i++){
		elem = abs(inp_vector[i]);
		if (elem > norm){
			norm = elem;
		}
	}
	return norm;
}

void LinearAlgebra::LU_factorization(){
    double sum, deno;
    for(int k=0; k<N1; k++){
        for (int j=k; j<N1; j++){
            sum = 0;
            for (int s=0; s<=k-1; s++){
                sum += Matrix[k*N2+s]*Matrix[s*N2+j]; //es lo mismo que L[k,s]*U[s,j]
            }
            Matrix[k*N2+j] -= sum; // es lo mismo que U[k,j]=Matrix[k,j]-sum
        }
        deno = Matrix[k*N2+k];
        if (deno==0.0){ //U[k,k]==0
            if(!pivot(k)){
                cout<<"No existe factorizacion LU"<<endl;
                return;
            }
            deno = Matrix[k*N2+k];
        }
        for (int i=k+1; i<N1; i++){
            sum =0;
            for (int s=0; s<=k-1; s++){
                sum+= Matrix[i*N2+s]*Matrix[s*N2+k]; //sum = L[i,k]*U[s,k]
            }
            Matrix[i*N2+k] = (Matrix[i*N2+k]- sum)/deno; //L[i,k] = (Matrix[i,k]-sum)/U[k,k]
        }
    }
    solve_method = 2;
}

void LinearAlgebra::LDL_factorization(){
    double sum;
    // Elementos de la diagonal Matriz[i,i] -> D[i,i]
    for(int j=0; j<N1; j++){ 
        sum = 0.0;
        for (int s=0; s<=j-1; s++){
            sum += Matrix[s*N2+s] * Matrix[j*N2+s] * Matrix[j*N2+s]; // sum += D[s,s]*L[j,s]^2
        }
        Matrix[j*N2+j] -= sum; //D[j,j] = Matrix[j,j]-sum
        if (Matrix[j*N2+j]==0){ //U[k,k]==0
            cout<<"No existe factorizacion LDL"<<endl;
            return;
        }
    //Elementos fuera de la diagonal, matriz simetrica Matrix[i,j]== Matrix[j,i] -> L[i,j]==L[j,i]
        for(int i = j+1; i<N1; i++){
            sum = 0.0;
            for (int s=0; s<=j-1; s++){
                sum += Matrix[s*N2+s] * Matrix[i*N2+s] * Matrix[j*N2+s]; //D[s,s]*L[i,s]*L[j,s]
            }
            Matrix[i*N2+j] = (Matrix[i*N2+j]-sum)/Matrix[j*N2+j]; //L[i,j] = (Matrix[i,j]-sum)/D[j,j]
            Matrix[j*N2+i] = Matrix[i*N2+j]; //L=L^T
        }
    }
    solve_method = 3;
}

void LinearAlgebra::solve(){
    double* new_vec = new double[N3]; //para conservar inp_intacto, 
    //se debe de crear un vector auxiliar que guarde su informacion
    //hasta que la solucion del sistema quede hecha en new_vec. Despues
    //new_vector será la direccion de memorio de out_vector miestras que el viejo será eliminado 
    switch (solve_method){
        case 1: //Eliminacion Gausiana
        elim_gausiana();
        break;

        case 2: //Factorizacion LU
        solve_triang_inf(1);
        change_ptr(&inp_vector, &out_vector); //inp->out out->inp //("la direccion de * esta en la variable *")
        change_ptr(&new_vec, &out_vector); //inp->new, new->out**
        solve_triang_sup();
        change_ptr(&new_vec,&inp_vector); //inp->inp** out->new**
        break;

        case 3: //Factorizacion LDL
        solve_triang_inf(1);
        change_ptr(&inp_vector, &out_vector); //inp->out out->inp ("la direccion de * esta en la variable *")
        change_ptr(&new_vec, &out_vector); //inp->new, new->out
        solve_diag();
        change_ptr(&inp_vector, &out_vector); //new->inp out->out**
        solve_triang_sup(1);
        change_ptr(&new_vec, &inp_vector); //inp->inp** new->new**
        break;
        
        case 4:
        solve_tridiagonal_system();
        break;

        case 5:  //Gradiente Conjugado
        conjugated_gradient(1000, 1e-8);
        break;

        case 6: //Simple multiplicacion de la matriz inversa
        solve_for_inv_method();
        break;
    }
    delete[] new_vec; //adios, gracias por ayudar :)
}

void LinearAlgebra::sum_rows(double* ext_vector){
    double sum;
    for (int i=0; i<N1; i++){
        sum=0;
        for(int j=0; j<N2; j++) sum+= Matrix[i*N2+j];
        ext_vector[i]=sum;
    }
}

void LinearAlgebra::sum_colms(double* ext_vector){
    double sum;
    for (int j=0; j<N2; j++){
        sum = 0;
        for(int i=0; i<N1; i++) sum+= Matrix[i*N2+j];
        ext_vector[j]=sum;
    }    
}

bool LinearAlgebra::is_diag_dom(){
    if(cache["diag_dom"] == 1) return false;
    if(cache["diag_dom"] == 2) return true;
    double sum;
    for (int i=0; i<N1; i++){
        sum=0;
        for(int j=0; j<N2; j++) sum+= Matrix[i*N2+j];
        if (abs(Matrix[i*N2+i])<abs(sum-Matrix[i*N2+i])){ 
            cache["diag_dom"] = 1;
            return false;
        }
    }
    cache["diag_dom"] = 2;
    return true;
}

void LinearAlgebra::jacobbi(int max_iter = 1000, double epsilon = 0.001){
    double sum, diag, error;
    double* temp = new double[N3];
    
    if (!is_diag_dom()) cout<<"Matriz no es diagonalmente dominante, es posible que el metodo de jacobbi no converja"<<endl;
    sum_rows(out_vector); //primer vector de iteracion, igual a la suma de las filas de la matrix dada
    for (int iter = 0; iter < max_iter; iter++){
        for (int i =0; i<N1; i++){
            sum = 0;
            diag = Matrix[i*N2+i];
            if (diag==0){
                cout<<"Metodo jacobbi no resoluble. Cero en diagonal"<<endl;
                delete[] temp;
                return;
            }
            for (int j=0; j<N2; j++) sum+=Matrix[i*N2+j]*out_vector[j];
            sum-= diag*out_vector[i];
            temp[i] = (inp_vector[i] -sum)/diag;
        }
        error = cuadratic_error(out_vector,N3, temp,N3);
        //intercambio de punteros
        change_ptr(&temp, &out_vector);

        if(error < epsilon){
            delete[] temp;
            return;
        }
    }
    cout<<"Metodo de jacobbi, no convergio en "<<max_iter<<" iteraciones"<<endl;
}

void LinearAlgebra::gauss_seidel(int max_iter = 1000, double epsilon = 0.001){
    double sum1, sum2, diag, error;
    double* temp = new double[N3];

    if (!is_diag_dom()) cout<<"Matriz no es diagonalmente dominante, es posible que el metodo gauss seidel no converja"<<endl;
    sum_rows(out_vector); //primer vector de iteracion, igual a la suma de las filas de la matrix dada
    for (int iter = 0; iter < max_iter; iter++){
        for (int i =0; i<N1; i++){
            sum1 = 0; 
            sum2 = 0;
            diag = Matrix[i*N2+i];
            if (diag==0){
                cout<<"Metodo gauss-seidel no resoluble. Cero en diagonal"<<endl;
                delete[] temp;
                return;
            }
            for (int j=0; j<i; j++) sum1+=Matrix[i*N2+j]*temp[j];
            for (int j=i+1; j<N2; j++) sum2+= Matrix[i*N2+j]*out_vector[j];
            temp[i] = (inp_vector[i]-sum1-sum2)/diag;
        }
        error = cuadratic_error(out_vector,N3,temp,N3);
        //intercambio de punteros
        change_ptr(&temp, &out_vector);

        if(error < epsilon){
            delete[] temp;
            return;
        }
    }
    cout<<"Metodo de gauss-seidel, no convergio en "<<max_iter<<" iteraciones"<<endl;
    delete[] temp;
}

void LinearAlgebra::set_Matrix_as_tridiagonal(double* vector1, double* vector2, double* vector3){
    down_diagonal = vector1;
    diagonal = vector2;
    up_diagonal = vector3;
    solve_method = 4;
}

void LinearAlgebra::solve_tridiagonal_system(){
    // La matriz tiene que ser necesariamente simetrica para que el sistema funcione
    double* L = new double[N3];
    double* D = new double[N3];
    double* auxptr = nullptr; //puntero auxiliar

    //Factorizacon LDL
    L[0] = 0, D[0] = diagonal[0];
    for(int i=1; i<N3; i++){
        L[i] = up_diagonal[i-1]/D[i-1];
        D[i] = diagonal[i] - D[i-1] * pow(L[i],2);
    }
    //Resolviendo el sistema tridiagonal
    //Triangular Inferior
    out_vector[0]=0;
    for (int i=0; i<N3; i++) out_vector[i] = inp_vector[i] - L[i]*out_vector[i];
    auxptr = inp_vector, inp_vector = out_vector, out_vector = auxptr;
    //Diagonal
    for(int i=0; i<N3; i++) out_vector[i] = inp_vector[i]/D[i];
    auxptr = inp_vector, inp_vector = out_vector, out_vector = auxptr;
    //Triangular Superior
    out_vector[N3-1] = inp_vector[N3-1];
    for (int i=N3-2; i>=0; i--) out_vector[i] = inp_vector[i] - L[i+1]*out_vector[i+1];

    delete[] L;
    delete[] D;
}
//-----------
void LinearAlgebra::change_ptr(double** ptr1, double** ptr2){
    double* aux = nullptr;
    aux = *ptr1, *ptr1 = *ptr2, *ptr2 = aux;
    return;
}

double LinearAlgebra::dot(double* vector1, double* vector2, int size){
    double sum = 0;
    for(int i=0; i<size; i++) sum+= vector1[i]*vector2[i];
    return sum;
}

double LinearAlgebra::max(double* vector, int size){
    double max = 0;
    for (int i=0; i<size; i++) if(vector[i]>max) max = vector[i];
    return max;
}

double LinearAlgebra::max_abs(double* vector, int size){
    double max = 0;
    for (int i=0; i<size; i++) if(abs(vector[i])>max) max = abs(vector[i]);
    return max;
}

double LinearAlgebra::min(double* vector, int size){
    double min;
    min = numeric_limits<double>::max();
    for (int i=0; i<size; i++) if(vector[i] < min) min = vector[i]; 
    return min;
}

double LinearAlgebra::min_abs(double* vector, int size){
    double min;
    min = numeric_limits<double>::max();
    for (int i=0; i<size; i++) if(abs(vector[i]) < min) min = abs(vector[i]); 
    return min;
}

void LinearAlgebra::divide_for(double* vector, double scalar, int size){
    for(int i = 0; i<size; i++) vector[i]/=scalar;
}

void LinearAlgebra::norm_vec(double* vector, int size){
    divide_for(vector, sqrt(dot(vector,vector, size)),size);
}

double* LinearAlgebra::ones(int size){
    double* ones_vector= new double[size];
    for (int i=0; i < size; i++) ones_vector[i]=1;
    return ones_vector;
}

void LinearAlgebra::pow_method(double* temp1, int max_iter = 1000, double epsilon = 0.001){
    if (N1!=N2 || N1!= N3) {
        cout<<"Matriz no es cuadrada"<<endl;
        return;
    }
    double* aux_ptr = nullptr;
    double* temp2 = new double[N3]; //este se convertirá en la nueva allocacion de eigenvector
    double lambda1 = 0;
    double lambda2 = 0; //inicializar eigenvalores
    double sum;

    //Metodo de la potencia
    for (int iter=0; iter<max_iter; iter++){

        for(int i=0; i<N1; i++){
            sum = 0;
            for (int j=0; j<N2; j++) sum += Matrix[i*N2+j] * temp1[j];
            temp2[i]= sum;
        }
        //lambda2 = dot(temp2, temp2, N3)/dot(temp2, temp1, N3);
        lambda2 = temp2[0]/temp1[0];
        if (abs(lambda2-lambda1) < epsilon){
            eigenvalues.emplace_back(lambda2);
            norm_vec(temp2, N3);
            eigenvectors.emplace_back(temp2);
            num_eigens++;
            delete[] temp1;
            return;
        }
        else{
            divide_for(temp2, max_abs(temp2, N3), N3);
            aux_ptr = temp1, temp1 = temp2, temp2 = aux_ptr; //cambiamos los punteros
            lambda1 = lambda2;
        }
    }

    cout<<"Metodo de potencia no convergio. Se retornara eigenvector vacio"<<endl;
    eigenvalues.emplace_back(0);
    delete[] temp2;
    delete[] temp1;
    return;
}

void LinearAlgebra::inverse_pow_method(double* temp1, int max_iter = 1000, double epsilon = 0.001){
    if (N1!=N2 || N1 != N3) {
        cout<<"Matriz no es cuadrada"<<endl;
        return;
    }
    double* aux_ptr = nullptr;
    double* temp2 = new double[N3]; //este se convertirá en la nueva allocacion de eigenvector
    double lambda1 = 0;
    double lambda2 = 0; //inicializar eigenvalores
    set_inp_vector(temp1, N3); //lo inicializamos como vector lleno de unos. default case

    //Metodo de la potencia inversa
    for (int iter=0; iter<max_iter; iter++){
        solve(); //solucion en out_vector
        //lambda2 = dot(out_vector, inp_vector, N3)/dot(out_vector, out_vector, N3);
        lambda2 = inp_vector[0]/out_vector[0];
        if ( abs(lambda2-lambda1) < epsilon){
            eigenvalues.emplace_back(lambda2);
            for (int i=0; i<N3; i++) temp2[i] = out_vector[i]; //copia de la solucion. out_vector se puede usar para otras cosas
            norm_vec(temp2, N3);
            eigenvectors.emplace_back(temp2);
            num_eigens++;
            return;
        }
        else{
            divide_for(inp_vector, max_abs(inp_vector,N3), N3);
            aux_ptr = inp_vector, inp_vector = out_vector, out_vector = aux_ptr; //cambiamos los punteros
            lambda1 = lambda2;
        }
    }
    cout<<"Metodo de potencia inverso no convergio. Se retornará eigenvector vacio"<<endl;
    eigenvalues.emplace_back(0);
    delete[] temp2;
    delete[] temp1;
    return;
}

void LinearAlgebra::getall_eigens(string method ="pow_method", int maxiter = 1000){
    num_eigens = 0;
    double* temp1 = new double[N3];
    sum_rows(temp1);//inicializar el eigenvector
    divide_for(temp1,max_abs(temp1,N3),N3); //para evitar numeros grandes
    if(method=="pow_method") pow_method(temp1);
    if(method=="inverse_pow_method") LU_factorization(), inverse_pow_method(temp1);
    //if(method=="jacobi") jacobi_eigens(maxiter);
}

void LinearAlgebra::minus_proy(double* vector1, double* vector2, int size){
    for(int i=0; i<size; i++) vector1[i]-= dot(vector1, vector2, size);
}
//--------------

void LinearAlgebra::jacobi_eigens(double** Q_eigen, double** P_eigen, bool save_eigens = true, int maxiter = 1000){
    double max, angle, ZERO = 1e-8;
    int index1, index2;
    double* Q_eigenvecs = new double[N1*N2];
    double* P_eigenvals = new double[N1*N2];

    //Copia de la matriz
    for(int i = 0; i <N1; i++){
        for (int j = 0; j<N2; j++){
            P_eigenvals[i*N2+j] = Matrix[i*N2+j];
            Q_eigenvecs[i*N2+j] = 0; //matriz unitaria
        }
        Q_eigenvecs[i*N2+i] = 1;
    }

    //METODO DE JACOBI PARA ENCONTRAR EIGENVALORES Y EIGENVECTORES
    for(int i=0; i<N1-1; i++){
        max = 0;
        for(int j = i+1; j<N2; j++){ // poner toda la fila i por debajo de ZERO
            if (max < abs(P_eigenvals[i*N2+j])){ //encontramos el maximo absoluto
                max = abs(P_eigenvals[i*N2+j]);
                index1 =  i, index2 = j; //indice fila, indice columna
            }
        }    

        if (max > ZERO){ //Rotacion
            rotation_jacobi_eigens(P_eigenvals,Q_eigenvecs,N1, N2, index1, index2);
        }
    }

    //Refinamiento de los eigenvectores y eigevalores (debido a que poner ceros en las otras filas puede quitar los ceros de otras filas)
    for(int iter = 0; iter < maxiter; iter++){ //damos un numero maxiter de intentos
        max = max_abs_Mat(P_eigenvals, N1, N2, index1, index2, true, false); //encuentra el maximo de la matriz simetrica y no consideres la diagonal
        if(max > ZERO){
            rotation_jacobi_eigens(P_eigenvals, Q_eigenvecs, N1, N2, index1, index2);
        }
        else { //Se ha conseguido la matriz diagonal con margen de error menor a ZERO
            //Gaurdamos los eigenvalores y eivectores en otra ubicacion accesible para el usuario
            if (save_eigens){
            double* aux;
            for(int j = 0; j < N2; j++){
                aux = new double[N1];
                for (int i = 0; i<N1; i++){
                    aux[j] =  Q_eigenvecs[i*N2+j];
                }
                eigenvectors.push_back(aux);
                eigenvalues.push_back(P_eigenvals[j*N2+j]);
            }
            num_eigens = N2;
            delete[] P_eigenvals;
            delete[] Q_eigenvecs;
            return;
            }

            else{
                *Q_eigen = Q_eigenvecs;
                *P_eigen = P_eigenvals;
                return;
            }
        }
    }
    cout<<"No se ha conseguido la diagonalizacion completa de la matriz con el metodo de Jacobi en "<<maxiter<<" iteraciones. Se retornara vector y eigenvalore vacios"<<endl;

    delete[] Q_eigenvecs;
    delete[] P_eigenvals;
}

void LinearAlgebra::rotation_jacobi_eigens(double* matrix, double* Rmatrix, int n1, int n2, int index1, int index2){
    double angle = 0.5*atan2( 2*matrix[index1*n2+index2], matrix[index1*n2+index1]-matrix[index2*n2+index2]);
    double* a = new double[n1];
    //Actualizacion de la matriz de eigenvectores (Producto matricial por la nueva matriz de rotación)
    
    //copia del vector
    for (int i =0; i<n1; i++) a[i] = Rmatrix[i*n2+index1];
    for (int k = 0; k <n2; k++){
        Rmatrix[k*n2+index1] = Rmatrix[k*n2+index1]*cos(angle)+Rmatrix[k*n2+index2]*sin(angle);
        Rmatrix[k*n2+index2] = Rmatrix[k*n2+index2]*cos(angle)-a[k]*sin(angle);
    }
    
    //Primera rotacion
    for (int i =0; i<n1; i++) a[i] = matrix[i*n2+index1];
    for (int k = 0; k <n1; k++){
        matrix[k*n2+index1] = matrix[k*n2+index1]*cos(angle)+matrix[k*n2+index2]*sin(angle);
        matrix[k*n2+index2] = matrix[k*n2+index2]*cos(angle)-a[k]*sin(angle);
    }

    //Segunda rotacion
    for (int i =0; i<n2; i++) a[i] = matrix[index1*n2+i];
    for(int k = 0; k<n2; k++){
        matrix[index1*n2+k] = matrix[index1*n2+k]*cos(angle)+matrix[index2*n2+k]*sin(angle);
        matrix[index2*n2+k] = matrix[index2*n2+k]*cos(angle)-a[k]*sin(angle);
    }
    delete[] a;
}

double LinearAlgebra::max_abs_Mat(const double *matrix, int n1, int n2, int& index1, int &index2, bool simetric =  false, bool diag = true){
    double max = 0;
    if (simetric && diag){ //si es simetrica e inclusimos la diagonal
        for(int i=0; i<n1-1; i++){
            for (int j=i; j<n2; j++){
                if (max < abs(matrix[i*n2+j])){
                    max = abs(matrix[i*n2+j]);
                    index1 = i, index2 = j;
                }
            }
        }
    }
    else if (simetric && !diag){ //simetrica y no incluimos la diagonal
        for(int i=0; i<n1-1; i++){
            for (int j=i+1; j<n2; j++){
                if (max < abs(matrix[i*n2+j])){
                    max = abs(matrix[i*n2+j]);
                    index1 = i, index2 = j;
                }
            }
        }        
    }
    else{ //Busqueda completa
        for(int i=0; i<n1-1; i++){
            for (int j = 0; j<n2; j++){
                if (max < abs(matrix[i*n2+j])){
                    max = abs(matrix[i*n2+j]);
                    index1 = i, index2 = j;
                }
            }
        }
    }
    return max;
}

// Tarea 9------------------
void LinearAlgebra::QR_factorization(){
    Q = new double[N1*N2];
    R = new double[N1*N2];
    double* aux1 = new double[N1];
    double* aux2 = new double[N1];
    double sum;
    
    //Copiando la matriz en Q y aux, R es la matriz cero
    for(int j = 0; j < N2; j++){
        for(int i = 0; i< N1; i++){
            Q[i*N2+j] = Matrix[i*N2+j];
            aux1[i] = Matrix[i*N2+j];
            aux2[i] =  Matrix[i*N2+j];
            R[i*N2+j] = 0;
        }

        //Normalizando la columna j de Q
        sum = 0;
        for(int i=0; i<N1; i++) sum += pow(Q[i*N2+j],2.0);
        sum = sqrt(sum);
        for(int i=0; i<N1; i++) Q[i*N2+j] /= sum; // esto es qj columna
 
        //Obteniedo los valores de R fuera de la diagonal
        for(int k=0; k<j; k++){
            sum = 0;
            for(int i=0; i<N1; i++) sum += Q[i*N2+j]*aux1[i]; //esto es qk^T aj
            R[k*N2+j] = sum; //esto es rkj
            sum = 0; 
            for (int i=0; i<N1; i++) aux2[i] -= R[k*N2+j]*Q[i*N2+k]; 
        }

        //Termino de la diagonal en R
        sum = 0;
        for (int i =0; i<N1; i++) sum+= pow(aux2[i], 2.0);
        R[j*N2+j] = sqrt(sum);
    }
    
    for (int i=0; i<N1; i++){
        cout<<"[";
        for (int j=0; j<N2; j++){
            cout<<Q[i*N2+j]<<", ";
        }
        cout<<"]"<<endl;
    }
    cout<<endl;

    for (int i=0; i<N1; i++){
        cout<<"[";
        for (int j=0; j<N2; j++){
            cout<<R[i*N2+j]<<", ";
        }
        cout<<"]"<<endl;
    }
    cout<<endl;
    delete[] aux1;
    delete[] aux2;
    solve_method = 4; //resolver sistema de ecuacione en modo QR
}

void LinearAlgebra::subspace_eigens(int DIM, int maxiter = 1000, double epsilon = 1e-8){
    if (N1 != N2){
        printf("No es posible hacer iteracion en el subespacio. Matriz no cuadrada de dimension (%d, %d) fue dada\n", N1,N2);
        return;
    }
    //print_array();
    double* sub_Matrix1 = new double[N1*DIM];
    double* sub_Matrix2 = new double[N1*DIM];
    double *B1 = new double[N1*DIM];
    double *B2 = new double[DIM*DIM];
    double *Q_eigenvecs, *P_eigenvals;
    int iter = 0, N;

    //Inicializadno matriz en el subespacio
    for(int i =0; i<N1; i++){
        for (int j = 0; j<DIM; j++) {
            sub_Matrix1[i*DIM+j] = 0;
            sub_Matrix2[i*DIM+j] = 0;
        }
    }
    for(int i =0; i<DIM; i++){
        sub_Matrix1[i*DIM+i] = 1;
        sub_Matrix2[i*DIM+i] = 1;
    }
    
    //METODO DE ITERACION EN EL SUBESPACIO CON JACOBI
    N = N1;
    N1 = DIM, N2 = DIM; //necesario para que funcione jacobi
    while(iter< maxiter){

        B1 = dot_mat(Matrix, N, N, sub_Matrix1, N, DIM); //La salida es B1 con dimensiones N1xDIM
        sub_Matrix2 = transpose(sub_Matrix1, N, DIM); //trasposicion de la matrix de subespacio
        B2 = dot_mat(sub_Matrix2, DIM, N, B1, N, DIM); //La salida es B2 con dimensiones DIMxDIM
        change_ptr(&B2,&Matrix);
        jacobi_eigens(&Q_eigenvecs, &P_eigenvals, false);
        change_ptr(&B2, &Matrix);
        //Resta matricial B1 = B2 - P_eigenvals 
        for (int i =0; i<DIM; i++){
            for(int j=0; j<DIM; j++) B1[i*DIM+j] = B2[i*DIM+j] - P_eigenvals[i*DIM+j];
        }
        
        if(frobenius_norm(B1, DIM, DIM, 2.0) < epsilon){ //se cumple el criterio

            //Guardamos los eigenvalores y los eigenvectores
            double* aux;
            for(int j = 0; j < DIM; j++){
                aux = new double[DIM];
                for (int i = 0; i < DIM; i++){
                    aux[j] =  Q_eigenvecs[i*DIM+j];
                }
                eigenvectors.push_back(aux);
                eigenvalues.push_back(P_eigenvals[j*DIM+j]);
            }
            num_eigens = DIM;
            N1 = N, N2 = N;
            delete[] P_eigenvals;
            delete[] Q_eigenvecs;
            delete[] sub_Matrix1;
            delete[] sub_Matrix2;
            delete[] B1;
            delete[] B2;
            return;
        }

        else{
            delete[] sub_Matrix2;
            sub_Matrix2 = dot_mat(sub_Matrix1, N, DIM, Q_eigenvecs, DIM, DIM);
            gram_schmit(sub_Matrix2, N, DIM);
            delete[] sub_Matrix1;
            sub_Matrix1 = dot_mat(Matrix, N, N, sub_Matrix2, N, DIM);
            gram_schmit(sub_Matrix1, N, DIM);
            iter++;
        }
    }

    cout<<"Metodo de iteracion en el subespacio no convergio en "<<maxiter<<" iteraciones"<<endl;
    N1 = N, N2 = N;
    delete[] P_eigenvals;
    delete[] Q_eigenvecs;
    delete[] sub_Matrix1;
    delete[] sub_Matrix2;
    delete[] B1;
    delete[] B2;
    return;
    
}

double* LinearAlgebra::dot_mat( const double* mat1, int sizef1, int sizec1, const double* mat2, int sizef2, int sizec2){
    if (sizec1 != sizef2){
        printf("dimensiones de matrices no concuerdan. (%d, %d) y (%d, %d) fueron dadas\n", sizef1, sizec2, sizef2, sizec2);
        return nullptr;
    }
    double* mat3 = new double[sizef1*sizec2];
    double sum;
    for (int i=0; i < sizef1; i++){
        for (int j=0; j < sizec2; j++){
            sum=0;
            for (int m = 0; m < sizec1; m++) sum +=  mat1[i*sizec1+m] * mat2[m*sizec2+j];
            mat3[i*sizec2 + j] = sum;
        }
    }
    return mat3;
}

double* LinearAlgebra::transpose(double* mat, int sizef, int sizec){
    double* aux = new double[sizef*sizec];
    for (int i = 0; i < sizef; i++){
        for (int j = 0; j<sizec; j++) aux[j*sizef+i] = mat[i*sizec+j];
    }
    return aux;
}

double LinearAlgebra::frobenius_norm(double* mat, int sizef, int sizec, double p = 2.0){
    double sum=0;
    for (int i =0; i<sizef; i++){
        for (int j=0; j<sizec; j++) sum += pow(mat[i*sizec+j], p);
    }
    return pow(sum, 1.0/p);
}

void LinearAlgebra::gram_schmit(double* mat, int sizef, int sizec){
    double* vec1 = new double[sizef];
    double* vec2 = new double[sizef];
    double sum, scalar;
    for (int j = 0; j < sizec; j++){
        for (int i=0; i < sizef; i++) vec1[i] = mat[i*sizec+j]; //Copia de la columna j

        //Resta de la proyeccion de los otros vectores
        for ( int k =0; k < j; k++){
            for (int i =0; i<sizef; i++) vec2[i] = mat[i*sizec+k]; //Copia de la columna k < j 
            scalar = dot(vec1, vec2, sizef); //producto punto de las columnas
            for (int i =0; i < sizef; i++) mat[i*sizec+j] -= scalar*vec2[i];
        }

        //Normalizacion de la columna j
        sum = 0;
        for (int i = 0; i<sizef; i++) sum += pow(mat[i*sizec+j],2.0);
        sum = sqrt(sum);
        for (int i =0; i< sizef; i++) mat[i*sizec+j] /=sum; 
    }
    delete[] vec1;
    delete[] vec2;
}

//TAREA 10--------------------------------------------------------------------------

void LinearAlgebra::conjugated_gradient(double maxiter = 5000, double epsilon = numeric_limits<double>::min()){
    if (N1 != N2){
        printf("Matriz no es cuadrada, no se puede resolver con el metodo de gradiente conjugado. Dimensiones (%d, %d) fueron dadas\n", N1, N2);
        return;
    }
    double *residual_vec,*conjugated_vec, *aux_vec;
    double delta, deltaOld, beta, alpha, error, k = 0;

    //METODO DEL GRADIENTE CONJUGADO
    for(int i = 0; i< N3; i++){
        out_vector[i] = 0;
    }

    residual_vec = dot_mat(Matrix, N1, N2, out_vector, N3, 1); //calcula el residuo
    for (int i = 0; i< N3; i++){
        residual_vec[i] = inp_vector[i] - residual_vec[i];
    }

    delta = sqrt(dot(residual_vec, residual_vec, N3));

    while(delta > epsilon && k < maxiter){
        if (k==0){
            conjugated_vec = copy(residual_vec, N3, 1);
        }
        else{
            beta = delta/deltaOld;
            for (int i = 0; i< N3; i++){
                conjugated_vec[i] = residual_vec[i] + beta*conjugated_vec[i];
            }
        }

        aux_vec = dot_mat(Matrix, N1, N2, conjugated_vec, N3, 1);
        alpha = delta/dot(conjugated_vec, aux_vec, N3);
        for(int i=0; i< N3; i++){
            out_vector[i] = out_vector[i] +alpha*conjugated_vec[i];
            residual_vec[i] = residual_vec[i] - alpha*aux_vec[i]; 
        }
        deltaOld = delta;
        delta = dot(residual_vec, residual_vec, N3);
        delete[] aux_vec;
        k++;
    }
    if(k == maxiter) {
        cout<<"Gradiente Conjugado no convergio en "<<maxiter<<" iteraciones"<<endl;
    }
    else{
        cout<<"Numero de iteraciones para convergencia: "<<k<<endl;
        error = sqrt(dot(residual_vec,residual_vec, N3))/N3;
        cout<<"El error residual es "<<error<<endl;
    }
    delete[] aux_vec;
    delete[] residual_vec;
    delete[] conjugated_vec;
    return;
}

double* LinearAlgebra::copy(double* mat, int sizef, int sizec){
    double* aux = new double[sizef*sizec];
    for (int i=0; i< sizef; i++){
        for (int j = 0; j<sizec; j++){
            aux[i*sizec+j] = mat[i*sizec+j];
        }
    }
    return aux;
}

void LinearAlgebra::add_mat(double* mat1, int sizef1, int sizec1, double* mat2, int sizef2, int sizec2, double* mat3){
    if(sizef1 != sizef2 || sizec1 != sizec2){
        printf("No es posble haber la suma de matrices. Dimensiones (%d, %d) y (%d, %d) fueron dadas\n", sizef1, sizec1, sizef2, sizec2);
        return;
    }

    for(int i=0; i<sizef1; i++){
        for (int j =0; j<sizec1; j++){
            mat3[i*sizec1+j] = mat1[i*sizec1+j] + mat2[i*sizec1+j];
        }
    }
    return;
}

void LinearAlgebra::subtraction_mat(double* mat1, int sizef1, int sizec1, double* mat2, int sizef2, int sizec2, double* mat3){
    if(sizef1 != sizef2 || sizec1 != sizec2){
        printf("No es posble haber la suma de matrices. Dimensiones (%d, %d) y (%d, %d) fueron dadas\n", sizef1, sizec1, sizef2, sizec2);
        return;
    }

    for(int i=0; i<sizef1; i++){
        for (int j =0; j<sizec1; j++){
            mat3[i*sizec1+j] = mat1[i*sizec1+j] - mat2[i*sizec1+j];
        }
    }
    return;
}

void LinearAlgebra::set_solve_method(string new_method){
    if (new_method == "inv_method") solve_method = 6;
    else cout<<"Metodo no encontrado"<<endl;
}

void LinearAlgebra::solve_for_inv_method(){
    delete[] out_vector;
    out_vector = dot_mat(Matrix, N1, N2, inp_vector, N3,1);
    return;
}