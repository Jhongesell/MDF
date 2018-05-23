// Ejemplo de una funcion para resolver la ecuacion diferencial parcial en 1D
// Uxx=-Pi*Pi*cos(Pix)
// xi<=U<=xf
// U(xi)=vi y U(xf)=vf

public class fdm1d {
   private int Visualiza;

   // Constructor
   public fdm1d() {
      Visualiza = 1;
   }

   // Resuelve Ax=b usando el metodo Jacobi
   public void Jacobi(double[][] A, double[] x, double[] b, int n, int iter) {
      int i, j, m;
      double sum;
      double[] xt = new double [n];

      for (m = 0; m < iter; m ++) {
         for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = 0; j < n; j ++) {
               if ( i == j) continue;
               sum += A[i][j] * x[j];
            }
            if (A[i][i] == 0.0) return;
            xt[i] = (1.0 / A[i][i]) * (b[i] - sum);
         }
         for (i = 0; i < n; i++) x[i] = xt[i];
      }
      if (Visualiza != 0) {
         System.out.println(" ");
         System.out.println("Matriz");
         for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
               System.out.print(A[i][j] + " ");
            }
            System.out.println(" ");
         }
         System.out.println(" ");
         System.out.println("b");
         for (i = 0; i < n; i++) {
            System.out.print(b[i] + " ");
         }
         System.out.println(" ");
         System.out.println("x");
         for (i = 0; i < n; i++) {
            System.out.print(x[i] + " ");
         }
         System.out.println(" ");
      }
   }


   // Resuelve Ax=b usando el metodo Gauus-Siedel
   public void Gauss_Siedel(double[][] A, double[] x, double[] b, int n, int iter) {
      int i, j, m;
      double sum;

      for (m = 0; m < iter; m ++) {
         for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = 0; j < n; j ++) {
               if ( i == j) continue;
               sum += A[i][j] * x[j];
            }
            if (A[i][i] == 0.0) return;
            x[i] = (1.0 / A[i][i]) * (b[i] - sum);
         }
      }
      if (Visualiza != 0) {
         System.out.println(" ");
         System.out.println("Matriz");
         for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
               System.out.print(A[i][j] + " ");
            }
            System.out.println(" ");
         }
         System.out.println(" ");
         System.out.println("b");
         for (i = 0; i < n; i++) {
            System.out.print(b[i] + " ");
         }
         System.out.println(" ");
         System.out.println("x");
         for (i = 0; i < n; i++) {
            System.out.print(x[i] + " ");
         }
         System.out.println(" ");
      }
   }


   // Lado derecho de la ecuacion diferencial parcial
   public static double LadoDerecho(double x) {
      double pi = 3.1415926535897932384626433832;
      return -pi * pi * java.lang.Math.cos(pi * x);
   }


   // Funcion Principal ....
   public static void main(String[] args) {

      double xi = -1.0;               // Inicio del diminio
      double xf = 2.0;                // Fin del dominio
      double vi = -1.0;               // Valor en la frontera xi
      double vf = 1.0;                // Valor en la frontera xf
      int n = 11;                     // Particion
      int N = n - 2;                  // Nodos interiores
      double h = (xf - xi) / (n - 1); // Incremento en la malla

      double[][] A = new double[N][N];  // Matriz A
      double[] b = new double[N];       // Vector b
      double[] x = new double[N];       // Vector x

      double R = 1 / (h * h);
      double P = -2 / (h * h);
      double Q = 1 / (h * h);

      // Primer renglon de la matriz A y vector b
      A[0][0] = P;
      A[0][1] = Q;
      b[0] = LadoDerecho(xi) - vi * R;
      // Renglones intermedios de la matriz A y vector b
      for (int i = 1; i < N - 1; i++) {
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
      fdm1d ejem = new fdm1d();
      ejem.Gauss_Siedel(A, x, b, N, 1000);
      ejem.Jacobi(A, x, b, N, 1000);
   }

}
