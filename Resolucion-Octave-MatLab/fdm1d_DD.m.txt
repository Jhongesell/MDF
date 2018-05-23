%{
 Metodo de diferencias Finitas una Dimension para resolver el problema
 con condiones Dirichlet 
   p(x)Uxx + q(x)Ux + r(x)U = f(x), con solucion analitica s(x)


 asi, para el problema particular
  -Uxx + Vux = 0   en 0 < x < 1 con U(0) = 0 y U(1) = 1   

  Escribir en consola
   v=@(x) 100.0;  % velocidad, posibles valore 1,25,100, etc 
   p=@(x) -1.0;
   q=@(x)  v(x);
   r=@(x) 0;
   f=@(x) 0;
   s=@(x) (exp(v(x)*x)-1.0)/(exp(v(x))-1.0);
   [error,A,b,u,x,V] = fdm1d_DD(p,q,r,f,0,1,0,1,11,1,s,1);

   
 asi, para el problema particular
  Uxx  = -pi*pi*cos(pi*x) en 0 < x < 1 con U(0)=1 y U(1)=-1
   
  Escribir en consola
   p=@(x) 1;
   q=@(x) 0;
   r=@(x) 0;
   f=@(x) -pi*pi*cos(pi*x);
   s=@(x) cos(pi*x);
   [error,A,b,u,x,V] = fdm1d_DD(p,q,r,f,0,1,1,-1,40,1,s,1);
 
 
 
 asi­, para el problema particular
  Uxx  = -pi*pi*cos(pi*x) en -1 < x < 2 con U(-1)=-1 y U(2)=1
   
  Escribir en consola
   p=@(x) 1;
   q=@(x) 0;
   r=@(x) 0;
   f=@(x) -pi*pi*cos(pi*x);
   s=@(x) cos(pi*x);
   [error,A,b,u,x,V] = fdm1d_DD(p,q,r,f,-1,2,-1,1,11,1,s,1);
 
 
 
%}



%% [error,A,b,u,x,V] = fdm1d_DD(p,q,r,f,s,xi,xf,vi,vf,n,grf,s,ssw)
%{
 Definicion del problema
    a(x)Uxx + b(x)Ux + c(x)U = f(x)

 p=@(x) ; Coeficiente para Uxx
 q=@(x) ; Coeficiente para Uxx
 r=@(x) ; Coeficiente para Uxx
 f=@(x) ; Lado derecho de la ecuacion
 s=@(x) ; Solucion analitica

 a	Inicio del dominio
 b	Fin del dominio
 vi	Valor de la condicion inicial 
 vf	Valor de la condicion final
 n	Numero de nodos en la particion
 grf Muestra las graficas con (1), otro valor no las muesta
 sws Se proporciona solucion analitica o no 
%}
function [error,A,b,u,x,V] = fdm1d_DD(p,q,r,f,xi,xf,vi,vf,n,grf,s,sws)
  if n < 3 
     return
  end

  % Numero de incognitas del problema
  tm = n -2;

  % Declaracion de la matriz y vectores de trabajo
  %A = sparse(tm,tm);
  A = zeros(tm,tm);   % Matriz de carga
  b = zeros(tm,1);    % Vector de carga
  u = zeros(tm,1);    % Vector de solucion

  x = zeros(n,1);     % Vector de coordenadas de la particion
  V = zeros(n,1) ;    % Vector solucion
  
  h = (xf-xi)/(n-1);
  h1 = h*h;

  % Llenado de los puntos de la malla 
  for i = 1: n,
     x(i) = xi + (i-1)*h;
  end
   
  % Llenado de la matriz y vector
  b(1) = f(xi) - p(xi)*(vi/h1);
  A(1,2) = p(x(1))/h1 - q(x(1))/(2.0*h);
  A(1,1) = ((-2.0 * p(x(1))) / h1) + r(x(1));
  for i=2:tm-1,
     A(i,i-1) = p(x(i))/h1 - q(x(i))/(2.0*h);
     A(i,i) = ((-2.0 * p(x(i))) / h1) + r(x(i));
     A(i,i+1) = p(x(i))/h1 + q(x(i))/(2.0*h);      
     b(i) = f(x(i));
  end
  A(tm,tm-1) = p(x(tm))/h1 - q(x(tm))/(2.0*h);
  A(tm,tm) = ((-2.0 * p(x(tm))) / h1) + r(x(tm));
  b(tm) = f(x(tm+1))-p(x(tm+1))*(vf/h1);
 
  % Soluciona el sistema
  u = inv(A)*b;

  % Copia la solucion obtenida del sistema lineal al vector solucion
  V(1) = vi;
  for i=1:tm,
     V(i+1) = u(i);
  end
  V(n) = vf;



  % Encuentra el error en norma infinita usando una particion de PAR = 10
  error = 0;
  if sws == 1 
    par = 10;
    for i = 1: n-1,
       inc = (x(i+1)-x(i))/par;
       for j = 1:par+1,
          px = x(i)+inc*(j-1);
          e = abs(s(px)-l(px,x(i),V(i),x(i+1),V(i+1)));
          if e > error
             error = e;
          end
       end
    end
  end

  % Revisa si se graficara 
  if grf == 1
    if sws == 1
      % Calcula la solucion analitica en la particion de calculo
      ua = zeros(n,1);
      for i = 1: n,
          ua(i) = s(x(i));
      end
    end
    % Graficar la solucion numerica
    plot(x,V,'o'); 
    hold

    
    % Grafica la solucion analitica en una particion de tamano xPart
    if sws == 1
      xPart = 1000;
      h = (xf-xi)/(xPart-1);
      xx = zeros(xPart,1);
      xa = zeros(xPart,1);
      for i = 1: xPart,
         xx(i) = xi + (i-1)*h;    
         xa(i) = s(xx(i));
      end
      plot(xx,xa);
    
      % Grafica el error
      figure(2);
      plot(x,V-ua);
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Evalua el punto x en la recta dada por los puntos (x1,y1) y (x2,y2), se usa para el calculo
% de la norma infinito

      function y = l(x,x1,y1,x2,y2)
   y = y1+((y2-y1)/(x2-x1))*(x-x1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



