// Ejemplo de una ecuación diferencial parcial en 1D
// -Uxx+Ux=0
// 0<U<1
// U(0)=0 y U(1)=1


// Resuelve Ax=b usando el metodo Jacobi
function x=Jacobi(A,b,N)
x=zeros(N,1);
xt=zeros(N,1);
eps = 1e-5;
iter = 100;
for it = 1: iter
   for i = 1: N
      sum = 0.0;
      for j = 1: N
         if i == j then
            continue;
         end
         sum = sum + A(i,j) * xt(j);
      end
      if A(i,i) == 0.0 then
         exit(1);
      end
      xt(i) = (b(i) - sum) / A(i,i);
   end
   sw = 0;
   for i = 1: N
      if abs(x(i) - xt(i)) >= eps then
         sw = 1;
      end
      x(i)=xt(i);
   end
   if sw == 0 then
      break;
   end
end
endfunction


a=0;           // Inicio dominio
c=1;           // Fin dominio
M=50;          // Partición
N=M-2;         // Nodos interiores
h=(c-a)/(M-1); // Incremento en la malla
Y0=0;          // Condición inicial en el inicio del dominio
Y1=1;          // Condición inicial en el fin del dominio
A=zeros(N,N);  // Matriz A
b=zeros(N);    // Vector b

P=2/(h^2);
Q=-1/(h^2)+1/(2*h);
R=-1/(h^2)-1/(2*h);

// Primer renglon de la matriz A y vector b
A(1,1)=P;
A(1,2)=Q;
b(1)=-Y0*R;
// Renglones intermedios de la matriz A y vector b
for i=2:N-1 
  A(i,i-1)=R;
  A(i,i)=P;
  A(i,i+1)=Q;
end
// Relglon final de la matriz A y vector b
A(N,N-1)=R;
A(N,N)=P;
b(N)=-Y1*Q;

// Resuleve el sistema lineal Ax=b
x=Jacobi(A,b,N);

// Prepara la graficación
xx=zeros(M,1);
for i=1:M
    xx(i)=a+h*(i-1);
end
yy=zeros(M,1);   
yy(1)=Y0;         // Condición inicial
for i=1:N
  yy(i+1)=x(i);
end
yy(M)=Y1;         // Condición inicial
// Grafica la solución de la Ecuación Diferencial Parcial en 1D
plot2d(xx,yy)


