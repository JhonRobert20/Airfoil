%Jhon Robert Matamoros Vitonera
%Airfoil
clear all;
c=9.2; %chord 
s=num2str(4612);
NACA=s; 
d1=str2double(s(1)); %primer escalar
d2=str2double(s(2));%segundo escalar
d34=str2double(s(3:4)); %tercer escalar
m=d1/100;
p=d2/10;
t=d34/100;
x=linspace(0, c, 200);
yt =5*t*c*(.2969*(sqrt(x/c))+-.1260*(x/c)+-.3516*(x/c).^2+.2843*(x/c).^3+-.1015*(x/c).^4);
mn=230;
xy=c;
centred=(c/200)*(mn-200)/2;
for k = 1:length(x)
      if x(k) <= p*c
          yc(k)=m*(x(k)/p^2)*(2*p-(x(k)/c));
          dx(k)=(2*m)/p^2*(p-(x(k)/c));
      elseif x(k) > p*c
          yc(k)=m*((c-x(k))/(1-p)^2)*(1+(x(k)/c)-(2*p));
          dx(k)=((2*m)/(1-p)^2)*(p-(x(k)/c));
      end
      %u de arriba, l de abajo (xu,yu) ; (xl,yl)
      theta=atan(dx(k));
         
      %Caso 1
      xu(k)=x(k)-yt(k)*sin(theta)+centred;
      yu(k)=yc(k)+yt(k)*cos(theta);
      xl(k)=x(k)+yt(k)*sin(theta)+centred;
      yl(k)=yc(k)-yt(k)*cos(theta);
      
end

%Constantes
rho=1.225;
patm=101350;
Lx=xy;
Ly=xy;
n=mn;
m=mn;
dx=Lx/m;
dy=Ly/n;
Xi=0:1:m;
Yi=0:1:n;

%generamos matriz 
A(:,1)=0;
A(:,n)=0;
A(n,:)=0;
A(1,:)=0;

%BOUNDARY CONDITIONS
Yisus=200/9.2;
Vmag=77.8;

% valor medio
A(:,:)=0.5*(Vmag*Ly);

%definimos borde izquierdo y derecho
j=1;
while(j<=n+1)
    B=(j-1)*Ly/n;
    A(j,1)=B*Vmag;
    A(j,m+1)=B*Vmag;
    if(j==n+1)
        
        may = A(j,m+1);
    end
        
    j=j+1;
end

%definimos arriba y abajo
i=1; 
while(i<=m+1)
A(1,i)=0;
A(n+1,i)=may;
i=i+1;
end

%Definimos Airfoil '4612'
j=2;
i=2;
Yn=m/2+1;

while (i<201)
    l1(i-1)=int32(xl(i)*Yisus);%Recta inferiorr
    l2(i-1)=int32(xu(i)*Yisus);%Recta superiorr
    l3(i-1)=int32(yl(i)*Yisus)+Yn;%Altura inferior
    l4(i-1)=int32(yu(i)*Yisus)+Yn;%Altura superior
    A(l4(i-1),l2(i-1))=NaN;
    j=1;
    while((l4(i-1)-j)>l3(i-1))
        A(l4(i-1)-j,l2(i-1))=NaN;
        j=j+1;
    end
    A(l3(i-1),l1(i-1))=NaN;
    i=i+1;
end
   
i=2;
while (i<m+1)
    j=2;
    while(j<n+1)
    if isnan(A(j,i+1))&& isnan(A(j,i-1))
        A(j,i)=NaN;
    end
    if isnan(A(j+1,i-1))&&isnan(A(j-1,i-1))
        A(j,i)=NaN;
    end
    j=j+1;
    end
    i=i+1;
end
        
%Function
t=1;
Iterations=10000;
while t<=Iterations
 
 i=2;
    while(i<m+1)
    j=2;
    while(j<n+1)
        if isnan(A(j,i))==false && isnan(A(j+1,i))==false && isnan(A(j-1,i))==false && isnan(A(j,i+1))==false && isnan(A(j,i-1))==false
            A(j,i)=(((dx)^2)*(A(j,i+1)+A(j,i-1))+((dy)^2)*(A(j+1,i)+A(j-1,i)))/(2*(((dy)^2)+((dx)^2)));
        end
    j=j+1;
    end
    i=i+1;
    end
t=t+1;
end





%definicion matrices para campos de velocidad y presion
uu(m+1,:)=Vmag;
uu(:,1)=Vmag;
uu(:,m+1)=Vmag;
uu(1,:)=Vmag;%%UU

vv(:,n+1)=0;
vv(1,:)=0;
vv(n+1,:)=0;
vv(:,1)=0;%VV

P(m+1,:)=patm;
P(:,1)=patm;
P(:,m+1)=patm;
P(1,:)=patm;

cp(n+1,:)=0;
cp(:,m+1)=0;
%llenamos matriz con campos de velocidad 

i=1;
while i<m+2
    j=2;
    while j<m+1 
        if j==1
            if ((isnan(A(j+1,i))==false))
            uu(j,i)=(((A(j+1,i))-(A(j,i)))/dy);
            end
        end
            
        if (j>1&&j<m)
           if ((isnan(A(j+1,i))==false)&&(isnan(A(j-1,i))==false))
            uu(j,i)=(((A(j+1,i))-(A(j-1,i)))/dy*2); %Usando centred  difference 
           end
        end
        if j==m
            if ((isnan(A(j+1,i))==false)&&(isnan(A(j-1,i))==false))
            uu(j,i)=(((A(j,i))-(A(j-1,i)))/dy);
            end
        end
            
        %Esto simplemente para para todo i, menos i=m
        if i==1
            if ((isnan(A(j+1,i))==false))
            vv(j,i)= -(((A(j,i+1))-(A(j,i)))/dx);
            end
        end
        if (i>1&&i<m)
           if ((isnan(A(j,i+1))==false)&&(isnan(A(j,i-1))==false))
             vv(j,i)= -(((A(j,i+1))-(A(j,i-1)))/dx*2); %Usando centreddifference 
           end
        end
        if i==m
            if ((isnan(A(j,i+1))==false)&&(isnan(A(j,i-1))==false))
             vv(j,i)= -(((A(j,i))-(A(j,i-1)))/dx); %Usando centreddifference 
            end
        end
          if (1<i && i<m+1)  
            if ((isnan(A(j+1,i))==false)&&(isnan(A(j-1,i))==false)&&(isnan(A(j,i+1))==false)&&(isnan(A(j,i-1))==false))
            P(j,i)=patm+(rho/2)*(Vmag^2-(sqrt((vv(j,i))^2+ (uu(j,i))^2))^2);
            cp(j,i)=(-patm+P(j,i))/(0.5*((340)^2)*rho);
            end
            if ((isnan(A(j+1,i)))&&(isnan(A(j-1,i)))&&(isnan(A(j,i+1)))&&(isnan(A(j,i-1)))||P(j,i)==0)
              P(j,i)=NaN;
            end
          end
    j=j+1;
    end
i=i+1;
end
i=1;
%Corrección de errores
while i<m+2%No puede, al tener vv=0, eliminar
    uu(1,i)=0;
    uu(m+1,i)=0;
    i=i+1;
end

vv(121,211)=0;
arriba(200)=0;
abajo(200)=0;
i=15;

while i<217
%     arriba(i)=cp(j,i);
    j=106;
    while j<116
        
       if (cp(j-1,i)==0||cp(j,i+1)==0||cp(j+1,i)==0)
          arriba(i-14)=cp(j,i);
          
          j=j+11;
       else
       j=j+1;
       end
    end
    i=i+1;
end

i=15;
while i<217
    j=138;
    while j>116
        
       if (cp(j-1,i)==0||cp(j,i+1)==0||cp(j+1,i)==0)
          abajo(i)=cp(j,i);
          j=j-52;
       else
       j=j-1;
       end
    end
    i=i+1;
end
suma(202)=0;
i=1;
cont=0;
sumaepica=0;
while i<203
 suma(i)=abajo(i)-arriba(i);
 sumaepica=suma(i)+sumaepica;
 if suma(i)<0
     cont=cont+1;
 end
 i=i+1;
end
resul=202-cont;
%CL
cl=((sqrt((sumaepica)^2))*dx*(15)*(1/9.2));
%En cuanto al lift, se conseidera no realizado correctamente

%PLOTS

%Stream Function
figure(1);
contourf(Xi,Yi,A,40)
title('Stream Funtion');
colorbar;

%Campo velocidades
figure(2);
quiver(Xi(1:10:end),Yi(1:10:end),uu(1:10:end,1:10:end),vv(1:10:end,1:10:end))
title('Campo de Velocidades');

%Campo presion
figure(3);
contourf(Xi,Yi,P,10)
title('Campo de Presiones');
colorbar;
 