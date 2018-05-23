// Ejemplo de una ecuaciï¿½n diferencial parcial en 1D
// Uxx=-Pi*Pi*cos(Pix)
// 0<=U<=1
// U(0)=1 y Ux(0.5)=-Pi


function y=LadoDerecho(x)
   y=-%pi*%pi*cos(%pi*x);
endfunction

function y=SolucionAnalitica(x)
  y=cos(%pi*x);
endfunction


a=0;           // Inicio dominio
c=0.5;         // Fin dominio
M=40;          // Particiï¿½n
N=M-1;         // Nodos interiores
h=(c-a)/(M-1); // Incremento en la malla
Y0=1;          // Condiciï¿½n Dirchlet inicial en el inicio del dominio 
Y1=-%pi;       // Condiciï¿½n Neumann inicial en el fin del dominio
A=zeros(N,N);  // Matriz A
b=zeros(N,1);    // Vector b


R=1/(h^2);
P=-2/(h^2);
Q=1/(h^2);


// Primer renglon de la matriz A y vector b
A(1,1)=P;
A(1,2)=Q;
b(1)=LadoDerecho(a)-Y0*R; // Frontera dirichlet
// Renglones intermedios de la matriz A y vector b
for i=2:N-1 
  A(i,i-1)=R;
  A(i,i)=P;
  A(i,i+1)=Q;
  b(i)=LadoDerecho(a+h*i);
end
// Relglon final de la matriz A y vector b
A(N,N-1)=-1/(h^2);
A(N,N)=-1/(h^2);
b(N)=Y1/h;



// Resuleve el sistema lineal Ax=b
x=inv(A)*b;

// Prepara la graficaciï¿½n
xx=zeros(M,1);
zz=zeros(M,1);
for i=1:M
  xx(i)=a+h*(i-1);
  zz(i)=SolucionAnalitica(xx(i));
end

yy=zeros(M,1);   
yy(1)=Y0;         // Condiciï¿½n inicial
for i=1:N
  yy(i+1)=x(i);
end

// Grafica la soluciï¿½n de la Ecuaciï¿½n Diferencial Parcial en 1D
plot2d(xx,[yy,zz])

yy(N)
zz(N)
