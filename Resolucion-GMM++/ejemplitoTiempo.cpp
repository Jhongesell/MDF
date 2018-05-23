// Metodo de Diferencias Finitas una Dimension para resolver el problema
// Uxx  = -pi*pi*cos(pi*x) , 0 < U < 1, U(0)=-1 y U(1)=1
// con solucion analitica: cos(pix)
// Para conocer mas del metodo ver:
// http://mmc.geofisica.unam.mx/edp/Ejemplitos/EcuacionesDiferencialesParciales/FDM/Introducci%C3%B3n%20al%20M%C3%A9todo%20de%20Diferencias%20Finitas%20y%20su%20Implementaci%C3%B3n%20Computacional.pdf
//
// Autor: Antonio Carrillo Ledesma
// http://mmc.geofisica.unam.mx/acl



#include <gmm/gmm.h>
#include <time.h>
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

   // Tomar tiempo de ejecucion
   clock_t inicio, final; 

   int M=11;               // Particion
   int N=M-2;              // Nodos interiores
   double a=0;             // Inicio dominio
   double c=1;             // Fin dominio
   double h=(c-a)/(M-1);   // Incremento en la malla
   double Y0=1.0;          // Condicion inicial en el inicio del dominio
   double Y1=-1.0;         // Condicion inicial en el fin del dominio


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

   

   inicio=clock(); 
   // LU para matrices densa
   gmm::dense_matrix<double> AA(N, N);
   gmm::copy(A,AA);
   gmm::lu_solve(AA, x, b); 
   final=clock(); 
   std::cout << "LU"<< std::endl<< x << gmm::endl; 
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);



   gmm::identity_matrix PS;   // Optional scalar product for cg
   gmm::identity_matrix PR;   // Optional preconditioner
   gmm::iteration iter(1E-10);// Iteration object with the max residu
   size_t restart = 50;       // restart parameter for GMRES


   for(i = 0; i < N; i++) x[i] = 0.0;
   inicio=clock(); 
   gmm::cg(A, x, b, PS, PR, iter); // Conjugate gradient
   final=clock(); 
   std::cout << "CGM"<< std::endl<< x << std::endl;
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);

/*
   for(i = 0; i < N; i++) x[i] = 0.0;
   inicio=clock(); 
   gmm::bicgstab(A, x, b, PR, iter); // BICGSTAB BiConjugate Gradient Stabilized
   final=clock(); 
   std::cout << "BICGSTAB"<< std::endl<< x << std::endl;
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);
*/


   for(i = 0; i < N; i++) x[i] = 0.0;
   inicio=clock(); 
   gmm::gmres(A, x, b, PR, restart, iter); // GMRES generalized minimum residual
   final=clock(); 
   std::cout << "GMRES"<< std::endl<< x << std::endl;
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);

/*
   for(i = 0; i < N; i++) x[i] = 0.0;
   inicio=clock(); 
   gmm::qmr(A, x, b, PR, iter); // Quasi-Minimal Residual method.
   final=clock(); 
   std::cout << "Quasi-Minimal"<< std::endl<< x << std::endl;
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);
*/

   for(i = 0; i < N; i++) x[i] = 0.0;
   // computation of a preconditioner (ILUT)
   inicio=clock(); 
   gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > Pre(A, 50, 1e-10);
   gmm::gmres(A, x, b, Pre, restart, iter);  // execute the GMRES algorithm
   final=clock(); 
   std::cout << "GMRES preconditiones ILUT"<< std::endl<< x << std::endl;
   printf("\nEl Tiempo de calculo es: %f (seg)\n\n\n\n",(final - inicio)/(double)CLOCKS_PER_SEC);


   std::cout << "Error en el calculo"<< std::endl;
   std::cout << a << "	" << abs(Y0 - SA(a))  << gmm::endl; // Valor de U en la frontera izquierda
   for(i = 0; i < N; i++)
   {
      std::cout << (i + 1)*h << "	" << abs(x[i] - SA((a + (i + 1)*h))) << gmm::endl;
   }
   std::cout << c << "	" << abs(Y1 - SA(c)) << gmm::endl; // Valor de U en la frontera derecha

   return 0;
}

