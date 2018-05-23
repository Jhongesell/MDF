% Ejemplo de una funcion para resolver la ecuacion diferencial parcial en 1D
% Uxx=-Pi*Pi*cos(Pix)
% xi<=U<=xf
% U(xi)=vi y U(xf)=vf

% Ejecuta la funcion
% [A,b,x]=fdm1d(-1,2,-1,1,50);


%    function [A,b,x] = fdm1d(xi,xf,vi,vf,n)
% Donde
%    xi    Inicio del domio
%    xf    Fin del dominio
%    vi    Valor de la condicion de frontera en xi
%    vf    valor de la condicion de frontera final en xf
%    n     Tamano de la particion
function [A,b,x] = fdm1d(xi,xf,vi,vf,n)

    N=n-2;            % Nodos interiores
    h=(xf-xi)/(n-1);  % Incremento en la malla
    A=zeros(N,N);   % Matriz A
    b=zeros(N,1);     % Vector b
    
    
    R=1/(h^2);
    P=-2/(h^2);
    Q=1/(h^2);
    
    
    % Primer renglon de la matriz A y vector b
    A(1,1)=P;
    A(1,2)=Q;
    b(1)=LadoDerecho(xi)-vi*R;
    % Renglones intermedios de la matriz A y vector b
    for i=2:N-1 
      A(i,i-1)=R;
      A(i,i)=P;
      A(i,i+1)=Q;
      b(i)=LadoDerecho(xi+h*(i-1));
    end
    % Renglon final de la matriz A y vector b
    A(N,N-1)=R;
    A(N,N)=P;
    b(N)=LadoDerecho(xi+h*N)-vf*Q;
    
    % Resuleve el sistema lineal Ax=b
    x=inv(A)*b;
    
    % Prepara la graficacion
    xx=zeros(n,1);
    zz=zeros(n,1);
    for i=1:n
      xx(i)=xi+h*(i-1);
      zz(i)=SolucionAnalitica(xx(i));
    end
    
    yy=zeros(n,1);   
    yy(1)=vi;         % Condicion inicial
    for i=1:N
      yy(i+1)=x(i);
    end
    yy(n)=vf;        % Condicion inicial
    
    % Grafica la solucion de la Ecuacion Diferencial Parcial en 1D
    plot(xx,[yy,zz]);
endfunction


function y=LadoDerecho(x)
   y=-pi*pi*cos(pi*x);
endfunction

function y=SolucionAnalitica(x)
  y=cos(pi*x);
endfunction

% [A,b,x]=fdm1d(-1,2,-1,1,30);
