// Metodo de Diferencias Finitas una Dimension para resolver el problema
// Uxx  = -pi*pi*cos(pi*x) , 0 < U < 1, U(0)=-1 y U(1)=1
// con solucion analitica: cos(pix)
// Para conocer mas del metodo ver:
// http://mmc.geofisica.unam.mx/edp/Ejemplitos/EcuacionesDiferencialesParciales/FDM/Introducci%C3%B3n%20al%20M%C3%A9todo%20de%20Diferencias%20Finitas%20y%20su%20Implementaci%C3%B3n%20Computacional.pdf
//
// Autor: Antonio Carrillo Ledesma
// http://mmc.geofisica.unam.mx/acl



#include <gmm/gmm.h>
#include <math.h>

const double  pi = 3.141592653589793;

// Lado derecho
double LD(double x)
{
   return ( -pi * pi * cos(pi * x));
}
 
// Solucion analitica
double SA(double x)
{
   return (cos(pi * x));
}


// Ejemplito del Metodo de Diferencias Finitas una Dimension para resolver el problema
//    Uxx  = -pi*pi*cos(pi*x) , 0 < U < 1, U(0)=-1 y U(1)=1
// con solucion analitica: cos(pix)
// usando GMM++
int main(void)
{
   int M=11;               // Particion
   int N=M-2;              // Nodos interiores
   double a=0;             // Inicio dominio
   double c=1;             // Fin dominio
   double h=(c-a)/(M-1);   // Incremento en la malla
   double Y0=1.0;          // Condicion inicial en el inicio del dominio
   double Y1=-1.0;         // Condicion inicial en el fin del dominio


   // Matriz densa
   gmm::dense_matrix<double> AA(N, N);


   // Matriz dispersa
   gmm::row_matrix< gmm::rsvector<double> > A(N, N);

   // Vectores
   std::vector<double> x(N), b(N);


   int i;
   double P = -2 / (h * h);
   double Q = 1 / (h * h);
   double R = 1 / (h * h);



   A(0, 0) = P; // Primer renglon de la matriz A y vector b
   A(0, 1) = Q;
   b[0] = LD(a + h) - (Y0 / (h * h));

   for(i = 1; i < N - 1; i++) // Renglones intermedios de la matriz A y vector b
   {
      A(i, i - 1) = R;
      A(i, i) = P;
      A(i, i + 1) = Q;
      b[i] = LD(a + (i + 1) * h);
   }

   A(N - 1, N - 2) = R; // Relglon final de la matriz A y vector b
   A(N - 1, N - 1) = P;
   b[N - 1] = LD(a + (i + 1) * h) - (Y1 / (h * h));

   // Copia la matriz dispersa a la densa para usarla en LU
   gmm::copy(A,AA);

   // Visualiza la matriz y el vector
   std::cout << "Matriz A"<< AA << gmm::endl;
   std::cout << "Vector b"<< b << gmm::endl;


   // LU para matrices densa
   gmm::lu_solve(AA, x, b);
   std::cout << "LU"<< x << gmm::endl;



   gmm::identity_matrix PS;   // Optional scalar product for cg
   gmm::identity_matrix PR;   // Optional preconditioner
   gmm::iteration iter(10E-6);// Iteration object with the max residu
   size_t restart = 50;       // restart parameter for GMRES



   gmm::cg(A, x, b, PS, PR, iter); // Conjugate gradient
   std::cout << "CGM"<< x << std::endl;



   gmm::bicgstab(A, x, b, PR, iter); // BICGSTAB BiConjugate Gradient Stabilized
   std::cout << "BICGSTAB"<< x << std::endl;



   gmm::gmres(A, x, b, PR, restart, iter); // GMRES generalized minimum residual
   std::cout << "GMRES"<< x << std::endl;



   gmm::qmr(A, x, b, PR, iter); // Quasi-Minimal Residual method.
   std::cout << "Quasi-Minimal"<< x << std::endl;



   // Visualiza la solucion numerica
   std::cout << "Solucion Numerica"<< std::endl;
   std::cout << a << "  " << Y0 << gmm::endl; // Valor de U en la frontera izquierda
   for(i = 0; i < N; i++)
   {
      std::cout << (i + 1)*h << "       " << x[i] << gmm::endl;
   }
   std::cout << c << "  " << Y1 << gmm::endl; // Valor de U en la frontera derecha


   // Visualiza la solucion analitica
   std::cout << "Solucion Analitica"<< std::endl;
   std::cout << a << "  " << SA(a) << gmm::endl; // Valor de U en la frontera izquierda
   for(i = 0; i < N; i++)
   {
      std::cout << (i + 1)*h << "       " << SA((a + (i + 1)*h)) << gmm::endl;
   }
   std::cout << c << "  " << SA(c) << gmm::endl; // Valor de U en la frontera derecha


   // Visualiza el error en valor absoluto en cada nodo
   std::cout << "Error en el calculo"<< std::endl;
   std::cout << a << "	" << abs(Y0 - SA(a))  << gmm::endl; // Valor de U en la frontera izquierda
   for(i = 0; i < N; i++)
   {
      std::cout << (i + 1)*h << "	" << abs(x[i] - SA((a + (i + 1)*h))) << gmm::endl;
   }
   std::cout << c << "	" << abs(Y1 - SA(c)) << gmm::endl; // Valor de U en la frontera derecha

   return 0;
}
