#include <iostream>
#include <array>
#include <complex>

const size_t Nc = 3;
const size_t Noffd = 2*Nc*Nc-Nc;
const size_t Ndiag = 2*Nc;

const size_t Niter = 20;

// Same data layout as is in the chroma clover linop
struct Triang6 {
    double  diag[2*Nc];
    std::complex<double>  offd[Noffd];
};

// Computes A = B * C
// Ugly code...
void tri6_mult(Triang6& A, const Triang6& B, const Triang6& C){
    // Diagonal Terms:
    A.diag[0] = B.diag[0]  * C.diag[0]  
              + B.offd[0].real()  * C.offd[0].real()  + B.offd[0].imag()  * C.offd[0].imag()
              + B.offd[1].real()  * C.offd[1].real()  + B.offd[1].imag()  * C.offd[1].imag() 
              + B.offd[2].real()  * C.offd[2].real()  + B.offd[2].imag()  * C.offd[2].imag() 
              + B.offd[3].real()  * C.offd[3].real()  + B.offd[3].imag()  * C.offd[3].imag() 
              + B.offd[4].real()  * C.offd[4].real()  + B.offd[4].imag()  * C.offd[4].imag(); 

    A.diag[1] = B.diag[1]  * C.diag[1]  
              + B.offd[0].real()  * C.offd[0].real()  + B.offd[0].imag()  * C.offd[0].imag()
              + B.offd[5].real()  * C.offd[5].real()  + B.offd[5].imag()  * C.offd[5].imag() 
              + B.offd[6].real()  * C.offd[6].real()  + B.offd[6].imag()  * C.offd[6].imag() 
              + B.offd[7].real()  * C.offd[7].real()  + B.offd[7].imag()  * C.offd[7].imag() 
              + B.offd[8].real()  * C.offd[8].real()  + B.offd[8].imag()  * C.offd[8].imag(); 

    A.diag[2] = B.diag[2]  * C.diag[2]  
              + B.offd[1].real()  * C.offd[1].real()  + B.offd[1].imag()  * C.offd[1].imag()
              + B.offd[5].real()  * C.offd[5].real()  + B.offd[5].imag()  * C.offd[5].imag() 
              + B.offd[9].real()  * C.offd[9].real()  + B.offd[9].imag()  * C.offd[9].imag() 
              + B.offd[10].real() * C.offd[10].real() + B.offd[10].imag() * C.offd[10].imag() 
              + B.offd[11].real() * C.offd[11].real() + B.offd[11].imag() * C.offd[11].imag(); 

    A.diag[3] = B.diag[3]  * C.diag[3]  
              + B.offd[2].real()  * C.offd[2].real()  + B.offd[2].imag()  * C.offd[2].imag()
              + B.offd[6].real()  * C.offd[6].real()  + B.offd[6].imag()  * C.offd[6].imag() 
              + B.offd[9].real()  * C.offd[9].real()  + B.offd[9].imag()  * C.offd[9].imag() 
              + B.offd[12].real() * C.offd[12].real() + B.offd[12].imag() * C.offd[12].imag() 
              + B.offd[13].real() * C.offd[13].real() + B.offd[13].imag() * C.offd[13].imag(); 

    A.diag[4] = B.diag[4]  * C.diag[4]  
              + B.offd[3].real()  * C.offd[3].real()  + B.offd[3].imag()  * C.offd[3].imag()
              + B.offd[7].real()  * C.offd[7].real()  + B.offd[7].imag()  * C.offd[7].imag() 
              + B.offd[10].real() * C.offd[10].real() + B.offd[10].imag() * C.offd[10].imag() 
              + B.offd[12].real() * C.offd[12].real() + B.offd[12].imag() * C.offd[12].imag() 
              + B.offd[14].real() * C.offd[14].real() + B.offd[14].imag() * C.offd[14].imag(); 

    A.diag[5] = B.diag[5]  * C.diag[5]  
              + B.offd[4].real()  * C.offd[4].real()  + B.offd[4].imag()  * C.offd[4].imag()
              + B.offd[8].real()  * C.offd[8].real()  + B.offd[8].imag()  * C.offd[8].imag() 
              + B.offd[11].real() * C.offd[11].real() + B.offd[11].imag() * C.offd[11].imag() 
              + B.offd[13].real() * C.offd[13].real() + B.offd[13].imag() * C.offd[13].imag() 
              + B.offd[14].real() * C.offd[14].real() + B.offd[14].imag() * C.offd[14].imag();                 

    // Off-diagonal Terms:
    A.offd[0] = B.diag[0]  * C.offd[0]  + B.offd[0]  * C.diag[1]
              + B.offd[1]  * std::conj(C.offd[5])  + B.offd[2]  * std::conj(C.offd[6])
              + B.offd[3]  * std::conj(C.offd[7])  + B.offd[4]  * std::conj(C.offd[8]); 

    A.offd[1] = B.diag[0]  * C.offd[1]  + B.offd[0]  * C.offd[5]
              + B.offd[1]  * C.diag[2]  + B.offd[2]  * std::conj(C.offd[9])
              + B.offd[3]  * std::conj(C.offd[10]) + B.offd[4]  * std::conj(C.offd[11]); 

    A.offd[2] = B.diag[0]  * C.offd[2]  + B.offd[0]  * C.offd[6]
              + B.offd[1]  * C.offd[9]  + B.offd[2]  * C.diag[3]
              + B.offd[3]  * std::conj(C.offd[12]) + B.offd[4]  * std::conj(C.offd[13]); 

    A.offd[3] = B.diag[0]  * C.offd[3]  + B.offd[0]  * C.offd[7]
              + B.offd[1]  * C.offd[10] + B.offd[2]  * C.offd[12]
              + B.offd[3]  * C.diag[4]  + B.offd[4]  * std::conj(C.offd[14]); 

    A.offd[4] = B.diag[0]  * C.offd[4]  + B.offd[0]  * C.offd[8]
              + B.offd[1]  * C.offd[11] + B.offd[2]  * C.offd[13]
              + B.offd[3]  * C.offd[14] + B.offd[4]  * C.diag[5];      




    A.offd[5] = std::conj(B.offd[0])  * C.offd[1]  + B.diag[1]  * C.offd[5]
              + B.offd[5]  * C.diag[2]  + B.offd[6]  * std::conj(C.offd[9])
              + B.offd[7]  * std::conj(C.offd[10]) + B.offd[8]  * std::conj(C.offd[11]);   

    A.offd[6] = std::conj(B.offd[0])  * C.offd[2]  + B.diag[1]  * C.offd[6]
              + B.offd[5]  * C.offd[9]  + B.offd[6]  * C.diag[3]
              + B.offd[7]  * std::conj(C.offd[12]) + B.offd[8]  * std::conj(C.offd[13]);   

    A.offd[7] = std::conj(B.offd[0])  * C.offd[3]  + B.diag[1]  * C.offd[7]
              + B.offd[5]  * C.offd[10] + B.offd[6]  * C.offd[12]
              + B.offd[7]  * C.diag[4]  + B.offd[8]  * std::conj(C.offd[14]);   

    A.offd[8] = std::conj(B.offd[0])  * C.offd[4]  + B.diag[1]  * C.offd[8]
              + B.offd[5]  * C.offd[11] + B.offd[6]  * C.offd[13]
              + B.offd[7]  * C.offd[14] + B.offd[8]  * C.diag[5];      




    A.offd[9]  = std::conj(B.offd[1])  * C.offd[2]  + std::conj(B.offd[5])  * C.offd[6]
               + B.diag[2]  * C.offd[9]  + B.offd[9]  * C.diag[3]
               + B.offd[10] * std::conj(C.offd[12]) + B.offd[11] * std::conj(C.offd[13]);   

    A.offd[10] = std::conj(B.offd[1])  * C.offd[3]  + std::conj(B.offd[5])  * C.offd[7]
               + B.diag[2]  * C.offd[10] + B.offd[9]  * C.offd[12]
               + B.offd[10] * C.diag[4]  + B.offd[11] * std::conj(C.offd[14]);   

    A.offd[11] = std::conj(B.offd[1])  * C.offd[4]  + std::conj(B.offd[5])  * C.offd[8]
               + B.diag[2]  * C.offd[11] + B.offd[9]  * C.offd[13]
               + B.offd[10] * C.offd[14] + B.offd[11] * C.diag[5];    




    A.offd[12] = std::conj(B.offd[2])  * C.offd[3]  + std::conj(B.offd[6])  * C.offd[7]
               + std::conj(B.offd[9])  * C.offd[10] + B.diag[3]  * C.offd[12]
               + B.offd[12] * C.diag[4]  + B.offd[13] * std::conj(C.offd[14]);   

    A.offd[13] = std::conj(B.offd[2])  * C.offd[4]  + std::conj(B.offd[6])  * C.offd[8]
               + std::conj(B.offd[9])  * C.offd[11] + B.diag[3]  * C.offd[13]
               + B.offd[12] * C.offd[14] + B.offd[13] * C.diag[5];   




    A.offd[14] = std::conj(B.offd[3])  * C.offd[4]  + std::conj(B.offd[7])  * C.offd[8]
               + std::conj(B.offd[10]) * C.offd[11] + std::conj(B.offd[12]) * C.offd[13]
               + B.diag[4]  * C.offd[14] + B.offd[14] * C.diag[5];    
}



// Computes Tr(A)
double tri6_trace(const Triang6& A){
    double tr = 0.0;
    for (size_t i = 0; i < Ndiag; i++)
        tr += A.diag[i];
    return tr;
}

// Computes Tr(A*B) using hermiticity of the matrices
double tri6_trace_mul(const Triang6& A, const Triang6& B){
    double trd, tro;
    for (size_t i = 0; i < Ndiag; i++)
        trd += A.diag[i] * B.diag[i];

    for (size_t i = 0; i < Noffd; i++)
        tro += A.offd[i].real() * B.offd[i].real() + A.offd[i].imag() * B.offd[i].imag();
    return trd + tro + tro;
}   

int main(){
    Triang6 A, A2, A3, tmp;
    int sign = 0;

    // Define A
    for (size_t i = 0; i < Ndiag; i++)
        A.diag[i] = 1.8 - i*0.1 - 1.55;    // this is just some test data...
    for (size_t i = 0; i < Noffd; i++)
        A.offd[i] = std::complex<double>{ 1.2 - i*0.2, 1.2 - i*0.2 - 0.1 }; 

    // Define A^2 and A^3
    tri6_mult(A2, A, A);
    tri6_mult(A3, A, A2);
    
    std::array<double, 7> tr;
    tr[0] = 0;
    tr[1] = 0;
    tr[2] = tri6_trace(A2);
    tr[3] = tri6_trace(A3);
    tr[4] = tri6_trace_mul(A2, A2);
    tr[5] = tri6_trace_mul(A2, A3);
    tr[6] = tri6_trace_mul(A3, A3);

    std::cout << "Tr = [ ";
    for (size_t i = 2; i < 7; i++)
        std::cout << tr[i] << "  ";
    std::cout << "] " << std::endl;
    
    
    std::array<double, 5> p;
    p[0] = -1.0/6.0*tr[6] + 1.0/18.0*tr[3]*tr[3] - 1.0/48.0*tr[2]*tr[2]*tr[2] + 1.0/8.0*tr[4]*tr[2];
    p[1] = -1.0/5.0*tr[5] + 1.0/6.0 *tr[2]*tr[3];
    p[2] = -1.0/4.0*tr[4] + 1.0/8.0 *tr[2]*tr[2];
    p[3] = -1.0/3.0*tr[3];
    p[4] = -1.0/2.0*tr[2];

    std::cout << "p = [ ";
    for (size_t i = 0; i < 5; i++)
        std::cout << p[i] << "  " ;
    std::cout << "] " << std::endl;

    std::array<double, Niter+1> c;
    c[0] = 1;
    for (size_t i = 1; i < Niter+1; i++)
        c[i] = c[i - 1] / i;
    if (sign == 1)
        for (size_t i = 1; i < Niter+1; i += 2)
            c[i] = -c[i];

    std::cout << "c = [ ";
    for (size_t i = 0; i < Niter+1; i++)
        std::cout << c[i] << "  " ;
    std::cout << "] " << std::endl;

    std::array<double, 6> q;    
    if (Niter > 5){ 
        int ic = Niter - 6;
        for (size_t i = 0; i < 6; i++)
            q[i] = c[ic+1+i];
        
        while (ic >= 0) {
            double q5 = q[5];
            q[5] = q[4];
            q[4] = q[3]  - q5 * p[4];
            q[3] = q[2]  - q5 * p[3];
            q[2] = q[1]  - q5 * p[2];
            q[1] = q[0]  - q5 * p[1];
            q[0] = c[ic] - q5 * p[0];
            ic -= 1;
        }
    }
    else {
        std::cout << "error" << std::endl;
    }
    std::cout << "q = [ ";
    for (size_t i = 0; i < 6; i++)
        std::cout << q[i] << "  ";
    std::cout << "] " << std::endl;

    // I*q0 + A*q1 + A^2*q2 + A^3*q3
    for (size_t i = 0; i < Ndiag; i++)
        A.diag[i] = q[0] + q[1]*A.diag[i] + q[2]*A2.diag[i] + q[3]*A3.diag[i];   
    for (size_t i = 0; i < Noffd; i++)
        A.offd[i] = q[1]*A.offd[i] + q[2]*A2.offd[i] + q[3]*A3.offd[i];

    // A^4*q4
    tri6_mult(tmp, A2, A2);
    for (size_t i = 0; i < Ndiag; i++)
        A.diag[i] += q[4]*tmp.diag[i];   
    for (size_t i = 0; i < Noffd; i++)
        A.offd[i] += q[4]*tmp.offd[i];

    // A^5*q5
    tri6_mult(tmp, A3, A2);
    for (size_t i = 0; i < Ndiag; i++)
        A.diag[i] += q[5]*tmp.diag[i];   
    for (size_t i = 0; i < Noffd; i++)
        A.offd[i] += q[5]*tmp.offd[i];


    int ii = 0;
    for (size_t i = 0; i < 6; i++){
        for (size_t j = 0; j < 6; j++){
            if(i == j) std::cout << A.diag[i] << "  \t\t";
            else if(i < j)  {std::cout << A.offd[ii] << "\t"; ii++; }
            else std::cout << "    " << 0 << " \t\t\t";
        }
        std:: cout << std::endl;    
    }
    return 0;
}
