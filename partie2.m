%%Question 1
load('projet_partie2.mat');
data = xlsread('quarterly.xls');
[n,k]=size(data);

%%Question 2
inflation = [];
CPI = quarterly(:,9);
for i = 1:211 
    inflation(i+1) = (CPI(i+1)-CPI(i))/CPI(i+1);
end
figure;
plot(inflation)

%%Question 3
autocorrelogramme = autocorr(inflation);

partiel = parcorr(inflation) ;

%%Question 4
%Cf rapport

%%Question 5
ar(inflation, 3);

ar1 = acf(1);
ar2 = acf(2);
ar3 = acf(3);
ar4 = acf(4);

inf_1 = inflation;
for (i = 4:212)
    inf_1(i) = ar1 + ar2 * inflation(i-1) + ar3 * inflation(i-2) + ar4 * inflation(i-3) ; 
end

figure
hold on;
plot ( 1:212 , inf_1) ;
title("courbe d'autoregression d'ordre 3");
hold off;

%%Question 6
[n,k]=size(data);
unemp = data(:,13);


phillips = plot(inflation, unemp, '.');
[fit2, gof2, fitinfo2] = fit(inflation.', unemp, 'a+b*x');



const = ones(n,1);

y = unemp;
X = [const,inflation];
[n,k]=size(X);
beta=inv(X'*X)*X'*y;
u=y-X*beta;

x = inf;
y = const*beta(1)+inf*beta(2);
figure;
hold on
plot(x,y)
scatter(inflation,unemp);
title('courbe de Philips');
hold off
xlabel('Inflation');
ylabel('Unemployement');


%%Question 7 : test de Durbin-Watson
ut_1 = u(1:n-1);
ut = u(2:n);
y=ut;
X=[ut_1];
[n,k]=size(X);
beta=inv(X'*X)*X'*y;
u=y-X*beta;
sig2=u'*u/(n-k);
std=sqrt(diag(sig2*inv(X'*X)));
t=(beta)./std;
fprintf('The test value of the autocorrelation is : %d \n',t);

pval = tdis_prb(t,n-k);
pval

%%Question 8
coefficients = fgls(inflation.', unemp, 'display', 'final');


%%Question 9
inflation_p1 = inflation(1:105);
inflation_p2 = inflation(106:211);

unemp_p1 = unemp(1:105);
unemp_p2 = unemp(106:211);

coeff_p1 = fgls(inflation_p1.', unemp_p1, 'display', 'final');
coeff_p2 = fgls(inflation_p2.', unemp_p2, 'display', 'final');

X_1 = [ones(105, 1), unemp_p1];
[n,k] = size(X_1);
beta = inv(X_1'*X_1)*X_1'*inflation_p1;
u = inflation_p1 - X_1*beta;
SSR0=u'*u;


X_2 = [ones(106, 1), unemp_p2];
[n,k] = size(X_2);
beta = inv(X_2.'*X_2)*X_2.'*inflation_p2;
v = inflation_p2 - X_2*beta;
SSR1=v'*v;


X = [ones(212, 1), unemp];
[n,k] = size(X_2);
beta = inv(X.'*X)*X.'*inflation;
w = unemp - X*beta;
SSR2=w'*w;


[n,k] = size(inf);
k = 2;
F = ((SSR0-(SSR1+SSR2))/k)/((SSR1+SSR2)*(n-2*k));
fprintf('The F value is : %d \n',F);

pval=fdis_prb(F,2,n-2*k);
fprintf('The p-value is : %d \n',pval);


%%Question 10 
[chow, pValue, stat, cValue] = chowtest(inflation.', unemp, 42);

%to find the breakpoint : first point for which the chowtest changes value
i=0;
while chow = 0
    i = = i +1;
    [chow, pValue, stat, cValue] = chowtest(inflation.', unemp, i);
    disp(i) ;
return i;

%Tracé de la courbe des maximums des F
[n,k] = size(inflation);
n_min = fix(0.15*n);
n_max = fix((1-0.15)*n);
F = [];
List = [];
for i = n_min:1:n_max+1 
   const = ones(n,1);
   X = [const];
   X_top = X(1:i);
   X_bot = X(i:212);
   Unemp_top = unemp(1:i) ; 
   Unemp_bot = unemp(i:212) ; 
   beta_top = inv(X_top'*X_top)*X_top'*Unemp_top;
   beta_bot = inv(X_bot'*X_bot)*X_bot'*Unemp_bot;
   u_bot = Unemp_bot - X_bot*beta_bot ;
   u_top = Unemp_top - X_top*beta_top ;
   SSR_top =u_top'*u_top ;
   SSR_bot =u_bot'*u_bot ;
   Fn = ((SSR0-(SSR_top+SSR_bot))/k)/((SSR_top+SSR_bot)/(n-2*2));
   F = [F Fn];
   List = [List i];
end
figure;
hold on
plot(List, F);
hold off
Top = List(find(F==max(F))) 

granger_cause(inflation, unemp, 0.05, 212);

%%Question 11
[n,k]=size(data);
unemp1 = unemp(1:n-5);
unemp2 = unemp(2:n-4);
unemp3 = unemp(3:n-3);
unemp4 = unemp(4:n-2);
unemp5 = unemp(5:n-1);

inf1 = inflation(1:n-5);
inf2 = inflation(2:n-4);
inf3 = inflation(3:n-3);
inf4 = inflation(4:n-2);
inf5 = inflation(5:n-1);

y = unemp1;
X = [inf2,inf3,inf4,inf5,unemp2,unemp3,unemp4,unemp5];
beta=inv(X'*X)*X'*y;
u=y-X*beta;
sig2=u'*u/(n-k);
std=sqrt(diag(sig2*inv(X'*X)));
t=(beta)./std;


%%Question 12
inf1 = inflation(1:n-5);
inf2 = inflation(2:n-4);
inf3 = inflation(3:n-3);
inf4 = inflation(4:n-2);
inf5 = inflation(5:n-1);

unemp1 = unemp(1:n-5);

%Cas 1 
y = unemp1;
X = [inf1];
beta=inv(X'*X)*X'*y;
One = sum(beta);
fprintf('The 1st value is : %d \n',One);

%Cas 2 
y = unemp1;
X = [inf1,inf2];
beta=inv(X'*X)*X'*y;
Two = sum(beta);
fprintf('The 2nd value is : %d \n',Two);

%Cas 3 
y = unemp1;
X = [inf1,inf2,inf3];
beta=inv(X'*X)*X'*y;
Three = sum(beta);
fprintf('The 3rd value is : %d \n',Three);

%Cas 4 
y = unemp1;
X = [inf1,inf2,inf3,inf4];
beta=inv(X'*X)*X'*y;
Four = sum(beta);
fprintf('The 4th value is : %d \n',Four);

%Cas 5 
y = unemp1;
X = [inf1,inf2,inf3,inf4,inf5];
beta=inv(X'*X)*X'*y;
Five = sum(beta);
fprintf('The 5th value is : %d \n',Five);

figure;
y = [One,Two,Three,Four,Five];
plot(y);

%%Question 13
unemp1 = unemp(1:n-5);
unemp2 = unemp(2:n-4);
unemp3 = unemp(3:n-3);
unemp4 = unemp(4:n-2);
unemp5 = unemp(5:n-1);

inf1 = inflation(1:n-5);
inf2 = inflation(2:n-4);
inf3 = inflation(3:n-3);
inf4 = inflation(4:n-2);
inf5 = inflation(5:n-1);

y = inf1;
X = [inf2,inf3,inf4,inf5,unemp2,unemp3,unemp4,unemp5];
beta=inv(X'*X)*X'*y;
u=y-X*beta;
sig2=u'*u/(n-k);
std=sqrt(diag(sig2*inv(X'*X)));
t=(beta)./std;
