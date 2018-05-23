// Ejemplo de una ecuaciï¿½n diferencial parcial en 1D
// Autor: Antonio Carrillo Ledesma
// -Uxx-k^2U=0
// 0<=U<=1
// U(0)=1 y Ux(1)=ikU(1)


TEST = 1; // (0) Diferencias finitas, (1) Diferencias finitas exactas segun Yau Shu Wong y Guangrui Li

function y=LadoDerecho(x)
   y=0.0;
endfunction

function y=SolucionAnalitica(x, k)
  //y=cos(k*x)+%i*sin(k*x);
  y=exp(%i*k*x);
endfunction


K = 150;
KK = K*K;

a=0;           // Inicio dominio
c=1;           // Fin dominio
M=300;         // Particiï¿½n
N=M-1;         // Nodos interiores
h=(c-a)/(M-1); // Incremento en la malla
Y0=1;          // Condiciï¿½n Dirchlet inicial en el inicio del dominio 
Y1=%i*K;          // Condiciï¿½n Neumann inicial en el fin del dominio
A=zeros(N,N);  // Matriz A
b=zeros(N);    // Vector b


if TEST = 0 then
   R=-1/(h^2);
   P=2/(h^2)-KK;
   Q=-1/(h^2);
else
   R=-1/(h^2);
   P=(2*cos(K*h)+(K*h)^2)/(h^2) - KK;
   Q=-1/(h^2);
end

// Primer renglon de la matriz A y vector b
A(1,1)=P;
A(1,2)=Q;
b(1)=LadoDerecho(a)-Y0*R; // Frontera dirichlet
// Renglones intermedios de la matriz A y vector b
for i=2:N-1 
  A(i,i-1)=R;
  A(i,i)=P;
  A(i,i+1)=Q;
  b(i)=LadoDerecho(a+h*(i-1));
end
// Relglon final de la matriz A y vector b
if TEST = 0 then
   A(N,N-1)=1/(h^2);
   A(N,N)=-1/(h^2)+ Y1/h;
   b(N)=LadoDerecho(c)/2;  
else  
   A(N,N-1)=1/(h^2);
   A(N,N)=-1/(h^2)+ %i*sin(K*h)/(h^2);
   b(N)=LadoDerecho(c)/2;
end



// Resuleve el sistema lineal Ax=b
x=inv(A)*b;

ESC = 5;
xxx=zeros(M*ESC,1);
zzz=zeros(M*ESC,1);
for i=1:M*ESC
  xxx(i)=a+h/ESC*(i-1);
  zzz(i)=SolucionAnalitica(xxx(i),K);
end

// Prepara la graficaciï¿½n
xx=zeros(M,1);
zz=zeros(M,1);
for i=1:M
  xx(i)=a+h*(i-1);
  zz(i)=SolucionAnalitica(xx(i),K);
end

yy=zeros(M,1);   
yy(1)=Y0;         // Condiciï¿½n inicial
for i=1:N
  yy(i+1)=x(i);
end

// Grafica la soluciï¿½n de la Ecuaciï¿½n Diferencial Parcial en 1D
plot2d(xx,yy,15)
plot2d(xxx,zzz)

