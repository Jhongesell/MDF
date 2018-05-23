// Ejemplo de una funcion para resolver la ecuacion diferencial parcial en 1D
// Uxx=-Pi*Pi*cos(Pix)
// xi<=U<=xf
// U(xi)=vi y U(xf)=vf


#include <math.h>
#include <stdio.h>

#define PARTICION 11  // Tamaño de la particion
#define N 9           // Numero de incognitas
#define VISUALIZA 1   // (0) No visualiza la salida, otro valor la visualiza


// Lado derecho de la ecuacion diferencial parcial
double LadoDerecho(double x)
{
   double pi = 3.1415926535897932384626433832;
   return -pi * pi * cos(pi * x);
}

// Resuelve Ax=b usando el metodo Jacobi
void Jacobi(double A[N][N], double x[], double b[], int n, int iter)
{
   int i, j, m;
   double sum;
   double xt[N];

   for (m = 0; m < iter; m ++)
   {
      for (i = 0; i < n; i++)
      {
         sum = 0.0;
         for (j = 0; j < n; j ++)
         {
            if ( i == j) continue;
            sum += A[i][j] * x[j];
         }
         if (A[i][i] == 0.0) return;
         xt[i] = (1.0 / A[i][i]) * (b[i] - sum);
      }
      for (i = 0; i < n; i++) x[i] = xt[i];
   }
   if (VISUALIZA)
   {
      printf("\nMatriz\n");
      for (i = 0; i < n; i++)
      {
         for (j = 0; j < n; j++)
         {
            printf("%f ", A[i][j]);
         }
         printf("\n");
      }
      printf("\nb\n");
      for (i = 0; i < n; i++)
      {
         printf("%f ", b[i]);
      }
      printf("\nx\n");
      for (i = 0; i < n; i++)
      {
         printf("%f ", x[i]);
      }
      printf("\n");
   }
}


// Resuelve Ax=b usando el metodo Gauus-Siedel
void Gauss_Siedel(double A[N][N], double x[], double b[], int n, int iter)
{
   int i, j, m;
   double sum;

   for (m = 0; m < iter; m ++)
   {
      for (i = 0; i < n; i++)
      {
         sum = 0.0;
         for (j = 0; j < n; j ++)
         {
            if ( i == j) continue;
            sum += A[i][j] * x[j];
         }
         if (A[i][i] == 0.0) return;
         x[i] = (1.0 / A[i][i]) * (b[i] - sum);
      }
   }
   if (VISUALIZA)
   {
      printf("\nMatriz\n");
      for (i = 0; i < n; i++)
      {
         for (j = 0; j < n; j++)
         {
            printf("%f ", A[i][j]);
         }
         printf("\n");
      }
      printf("\nb\n");
      for (i = 0; i < n; i++)
      {
         printf("%f ", b[i]);
      }
      printf("\nx\n");
      for (i = 0; i < n; i++)
      {
         printf("%f ", x[i]);
      }
      printf("\n");
   }
}



int main()
{

   double xi = -1.0;               // Inicio del diminio
   double xf = 2.0;                // Fin del dominio
   double vi = -1.0;               // Valor en la frontera xi
   double vf = 1.0;                // Valor en la frontera xf
   int n = PARTICION;              // Particion
   double h = (xf - xi) / (n - 1); // Incremento en la malla
   int i;

   double A[N][N];  // Matriz A
   double b[N];     // Vector b
   double x[N];     // Vector x

   double R = 1 / (h * h);
   double P = -2 / (h * h);
   double Q = 1 / (h * h);

   // Primer renglon de la matriz A y vector b
   A[0][0] = P;
   A[0][1] = Q;
   b[0] = LadoDerecho(xi) - vi * R;
   // Renglones intermedios de la matriz A y vector b
   for (i = 1; i < N - 1; i++)
   {
      A[i][i - 1] = R;
      A[i][i] = P;
      A[i][i + 1] = Q;
      b[i] = LadoDerecho(xi + h * i);
   }
   // Renglon final de la matriz A y vector b
   A[N - 1][N - 2] = R;
   A[N - 1][N - 1] = P;
   b[N - 1] = LadoDerecho(xi + h * N - 2) - vf * Q;

   // Resuleve el sistema lineal Ax=b
   Gauss_Siedel(A, x, b, N, 1000);
   Jacobi(A, x, b, N, 1000);
   return 0;
}

