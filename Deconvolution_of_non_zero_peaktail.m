function []=test2_with_fminsearch_tail()
clc;
close all;

Exp_matrix=load('20K_per_minute.txt.txt'); %%%% change the value for every beta
n=size(Exp_matrix,1);
c=0.1002;
Beta= [20/60];  %%%%% change the value for every beta

p=1;
for i=1:10:n   %%%%% change the middle value for every beta
    if (Exp_matrix(i,2)>=400) && (Exp_matrix(i,2)<=550)
    T_exp(p)=Exp_matrix(i,2)-c*Beta(1)+273.15;
    DVexpDT(p)=Exp_matrix(i,4)/100;
    p=p+1;
    end
end
p=p-1;

E_inf=90;
T=linspace(T_exp(1),T_exp(p),p);
E=linspace(0,E_inf,1.5*E_inf);    


k3=10^18.4445;
Lk3=log10(k3);

E3=65.8001;

sig3=7.415;

c3=0.1523;

b_guess = [Lk3 E3 sig3 c3 E_inf];
b_new=zeros(4);

options=optimset('PlotFcns',@optimplotfval,'MaxIter',100);
b_new = fminsearch(@(b)fun(b,Beta,T,DVexpDT), b_guess, options);

out=b_new.';
fid=fopen('output.txt','w');
fwrite(fid, b_new);
fclose(fid);
b_new

k3=10^b_new(1);

E3=b_new(2);

sig3=b_new(3);

c3=b_new(4);

%E_inf=b_new(1);

%E=linspace(0,E_inf,1.5*E_inf);

R=0.00198558;

pd3=makedist('normal',E3,sig3);

F3=pdf(pd3,E);


beta=Beta(1);    

for j=1:length(T)
    for i=1:length(E)
        f=(@(t)exp(-E(i)./R./(t)));
        g(j,i)=quadgk(f,1,T(j));
        y3(j,i)=(k3/beta)*exp(-E(i)/R/(T(j))-k3/beta*g(j,i))*F3(i);
    end
end

s=0;
for j=1:length(T)
    Y3(:)=y3(j,:);
    
    dadT3(j)=trapz(E,Y3);
    DVcalDT(j)=c3*dadT3(j);
    Res(j)=DVcalDT(j)-DVexpDT(j);
    s=s+(DVcalDT(j)-DVexpDT(j))^2;
end
plot(T,DVcalDT,T,DVexpDT,'k>');
end

function err=fun(b,Beta,T,DVexpDT)

k3=10^b(1);

E3=b(2);

sig3=b(3);

c3=b(4);

E_inf=b(5);

E=linspace(0,E_inf,1.5*E_inf);

R=0.00198558;

pd3=makedist('normal',E3,sig3);

F3=pdf(pd3,E);

beta=Beta(1);    

for j=1:length(T)
    for i=1:length(E)
        f=(@(t)exp(-E(i)./R./(t)));
        g(j,i)=quadgk(f,1,T(j));
        
        y3(j,i)=(k3/beta)*exp(-E(i)/R/(T(j))-k3/beta*g(j,i))*F3(i);
    end
end

s=0;
for j=1:length(T)
  
    Y3(:)=y3(j,:);

    dadT3(j)=trapz(E,Y3);
    DVcalDT(j)=c3*dadT3(j);
    s=s+(DVcalDT(j)-DVexpDT(j))^2;
end

err=s;
err=err(:);
%plot(T,DVcalDT,T,r1,'--',T,r2,'--',T,r3,'--');
%hold on;
%AllBetaData(:,kk)=DVcalDT(:);
%plot(T,AllBetaData);
end