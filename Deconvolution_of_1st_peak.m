function []=test2_with_fminsearch ()
clc;
close all;

Exp_matrix=load('20K_per_minute.txt');
n=size(Exp_matrix,1);
c=0.1002; % temp correction for the temperature lag, it will vary from TGA reactor to reactor.
Beta= [20/60];  

p=1;
for i=1:15:n  %%% 14=skipping points for 15C/min
    if (Exp_matrix(i,2)>=200) && (Exp_matrix(i,2)<=300)
    T_exp(p)=Exp_matrix(i,2)-c*Beta(1)+273.15;
    DVexpDT(p)=Exp_matrix(i,4)/100;
    p=p+1;
    end
end
p=p-1
    
T=linspace(T_exp(1),T_exp(p),p);
E=linspace(0,92,123);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the second value is E_inf, higher round value, The 3rd term should be 1.5 time of the second value
                                   
k1=10^15.9173;
Lk1=log10(k1);
k2=10^13.7779;
Lk2=log10(k2);
k3=10^18.8989;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modify from the tail opt%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lk3=log10(k3);

E1=40.7393;  
E2=45.5713;
E3=66.3500;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modify from the tail opt%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig1=0.0078;
sig2=1.4784;
sig3=6.9970;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modify from the tail opt%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1=0.0106;
c2=0.3688;
c3=0.1473;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modify from the tail opt%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_guess = [Lk1 Lk2 E1 E2 sig1 sig2 c1 c2]; %%%%%%%%%%%%%%initial values%%%%%%%%%%%%%%%%%%%
b_3 = [Lk3 E3 sig3 c3];   %%%%%%%%%%%%%%%%%%%%%%%constant%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=optimset('PlotFcns',@optimplotfval,'MaxIter',600);
b_new = fminsearch(@(b)fun(b,Beta,T,E,DVexpDT,b_3), b_guess, options);

out=b_new.';
fid=fopen('output.txt','w');
fwrite(fid, b_new);
fclose(fid);

BBBB=[b_new(2) b_new(4) b_new(6) b_new(8)] %%%% 

k1=10^b_new(1);
k2=10^b_new(2);


E1=b_new(3);
E2=b_new(4);


sig1=b_new(5);
sig2=b_new(6);


c1=b_new(7);
c2=b_new(8);


R=0.00198558;

pd1=makedist('normal',E1,sig1);
pd2=makedist('normal',E2,sig2);
pd3=makedist('normal',E3,sig3);

F1=pdf(pd1,E);
F2=pdf(pd2,E);
F3=pdf(pd3,E);


beta=Beta(1);    

for j=1:length(T)
    for i=1:length(E)
        f=(@(t)exp(-E(i)./R./(t)));
        g(j,i)=quadgk(f,1,T(j));
        
        y1(j,i)=(k1/beta)*exp(-E(i)/R/(T(j))-k1/beta*g(j,i))*F1(i);
        y2(j,i)=(k2/beta)*exp(-E(i)/R/(T(j))-k2/beta*g(j,i))*F2(i);
        y3(j,i)=(k3/beta)*exp(-E(i)/R/(T(j))-k3/beta*g(j,i))*F3(i);
    end
end

s=0;
for j=1:length(T)
    Y1(:)=y1(j,:);
    Y2(:)=y2(j,:);
    Y3(:)=y3(j,:);
    
    dadT1(j)=trapz(E,Y1);
    dadT2(j)=trapz(E,Y2);
    dadT3(j)=trapz(E,Y3);
    DVcalDT(j)=c1*dadT1(j)+c2*dadT2(j)+c3*dadT3(j);
    Res(j)=DVcalDT(j)-DVexpDT(j);
    s=s+(DVcalDT(j)-DVexpDT(j))^2;
end
plot(T,DVcalDT,T,DVexpDT,'k>');
end

function err=fun(b,Beta,T,E,DVexpDT,b_3)
k1=10^b(1);
k2=10^b(2);
k3=10^b_3(1);

E1=b(3);
E2=b(4);
E3=b_3(2);

sig1=b(5);
sig2=b(6);
sig3=b_3(3);

c1=b(7);
c2=b(8);
c3=b_3(4);

R=0.00198558;

pd1=makedist('normal',E1,sig1);
pd2=makedist('normal',E2,sig2);
pd3=makedist('normal',E3,sig3);

F1=pdf(pd1,E);
F2=pdf(pd2,E);
F3=pdf(pd3,E);


beta=Beta(1);    

for j=1:length(T)
    for i=1:length(E)
        f=(@(t)exp(-E(i)./R./(t)));
        g(j,i)=quadgk(f,1,T(j));
        
        y1(j,i)=(k1/beta)*exp(-E(i)/R/(T(j))-k1/beta*g(j,i))*F1(i);
        y2(j,i)=(k2/beta)*exp(-E(i)/R/(T(j))-k2/beta*g(j,i))*F2(i);
        y3(j,i)=(k3/beta)*exp(-E(i)/R/(T(j))-k3/beta*g(j,i))*F3(i);
    end
end

s=0;
for j=1:length(T)
    Y1(:)=y1(j,:);
    Y2(:)=y2(j,:);
    Y3(:)=y3(j,:);
    
    dadT1(j)=trapz(E,Y1);
    dadT2(j)=trapz(E,Y2);
    dadT3(j)=trapz(E,Y3);
    DVcalDT(j)=c1*dadT1(j)+c2*dadT2(j)+c3*dadT3(j);
    s=s+(DVcalDT(j)-DVexpDT(j))^2;
end

err=s;
err=err(:);
%plot(T,DVcalDT,T,r1,'--',T,r2,'--',T,r3,'--');
%hold on;
%AllBetaData(:,kk)=DVcalDT(:);
%plot(T,AllBetaData);
end