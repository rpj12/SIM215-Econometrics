% Partie 1
addpath(genpath('75jplv6'));

%% Question 1

load mroz.raw;
wage = mroz(:,7);
wagepositif = (wage>0);
mroz = mroz(wagepositif,:);


%% Question 2

age = mroz(:,5);
educ = mroz(:,6);
wage = mroz(:,7);

huswage = mroz(:,12);
medianhuswage = median(huswage);
huswagesup = (huswage>medianhuswage);
mrozsup = mroz(huswagesup,:);
huswageinf = (huswage<=medianhuswage);
mrozinf = mroz(huswageinf,:);

agesup = mrozsup(:,5);
ageinf = mrozinf(:,5);
wagesup = mrozsup(:,7);
wageinf = mrozinf(:,7);
educsup = mrozsup(:,6);
educinf = mrozinf(:,6);

mean(wage);
mean(wageinf);
mean(wagesup);
min(wage);
min(wageinf);
min(wagesup);
max(wage);
max(wageinf);
max(wagesup);
std(wage);
std(wageinf);
std(wagesup);
cov(wage);
cov(wageinf);
cov(wagesup);
corrcoef(huswage,wage);
corrcoef(huswage(huswageinf,:),wageinf);
corrcoef(huswage(huswagesup,:),wagesup);

mean(age);
mean(ageinf);
mean(agesup);
min(age);
min(ageinf);
min(agesup);
max(age);
max(ageinf);
max(agesup);
std(age);
std(ageinf);
std(agesup);
cov(age);
cov(ageinf);
cov(agesup);
corrcoef(huswage,age);
corrcoef(huswage(huswageinf,:),ageinf);
corrcoef(huswage(huswagesup,:),agesup);

mean(educ);
mean(educinf);
mean(educsup);
min(educ);
min(educinf);
min(educsup);
max(educ);
max(educinf);
max(educsup);
std(educ);
std(educinf);
std(educsup);
cov(educ);
cov(educinf);
cov(educsup);
corrcoef(huswage,educ);
corrcoef(huswage(huswageinf,:),educinf);
corrcoef(huswage(huswagesup,:),educsup);


%% Question 3

figure('Name','Histogramme de wage','NumberTitle','off');
histogram(wage);

lwage = mroz(:,21);
figure('Name','Histogramme de log(wage)','NumberTitle','off');
histogram(lwage);


%% Question 4

motheduc = mroz(:,15);
fatheduc = mroz(:,16);
disp('Corrélation entre motheduc et fatheduc');
disp(corrcoef(motheduc,fatheduc));


%% Question 5

figure('Name','Nuage de points de wage en fonction de educ','NumberTitle','off');
scatter(educ,wage,'filled');

figure('Name','Nuage de points de wage en fonction de exper','NumberTitle','off');
exper = mroz(:,19);
scatter(exper,wage,'filled');

figure('Name','Nuage de points de wage en fonction de fatheduc','NumberTitle','off');
scatter(fatheduc,wage,'filled');


%% Question 7

city = mroz(:,18);
nwifeinc = mroz(:,20);
kidslt6 = mroz(:,3);
kidsge6 = mroz(:,4);
[n,k] = size(mroz);
X7 = [ones(n,1) city educ exper nwifeinc kidslt6 kidsge6];
y7 = wage;
beta7 = inv(X7'*X7)*X7'*y7; 
u7 = y7-X7*beta7;
figure('Name','Histogramme des résidus','NumberTitle','off');
histogram(u7);

[n,k] = size(X7);
sig7 = u7'*u7/(n-k);
std7 = sqrt(diag(sig7*inv(X7'*X7)));
t7 = beta7./std7;

disp('coefficients de la régression linéaire de wage : ')
disp(beta7)


%% Question 8

[n,k] = size(mroz);
X8 = [ones(n,1) city educ exper nwifeinc kidslt6 kidsge6];
y8 = lwage;
beta8 = inv(X8'*X8)*X8'*y8;
u8 = y8-X8*beta8;
figure('Name','Histogramme des résidus','NumberTitle','off');
histogram(u8);

[n,k] = size(X8);
sig8 = u8'*u8/(n-k);
std8 = sqrt(diag(sig8*inv(X8'*X8)));
t8 = beta8./std8;

disp('coefficients de la régression linéaire de log(wage) : ')
disp(beta8)


%% Question 9

estim9 = t8(5);

if estim9<=tdis_inv(0.995,n-k) && estim9>=tdis_inv(0.005,n-k)
    disp('on accepte l hypothèse de non significativité de nwifeinc à 1%');
else
    disp('on rejette l hypothèse de non significativité de nwifeinc à 1%');
end

if estim9<=tdis_inv(0.975,n-k) && estim9>=tdis_inv(0.025,n-k)
    disp('on accepte l hypothèse de non significativité de nwifeinc à 5%');
else
    disp('on rejette l hypothèse de non significativité de nwifeinc à 5%');
end

if estim9<=tdis_inv(0.95,n-k) && estim9>=tdis_inv(0.05,n-k)
    disp('on accepte l hypothèse de non significativité de nwifeinc à 10%');
else
    disp('on rejette l hypothèse de non significativité de nwifeinc à 10%');
end

disp('p-value : ');
disp(tdis_prb(estim9,n-k));


%% Question 10

t10 = (beta8-0.01)./std8;
estim10 = t10(5);

if estim9<=tdis_inv(0.975,n-k) && estim9>=tdis_inv(0.025,n-k)
    disp('on accepte l hypothèse que coefnwifeinc = 0.01 à 5%');
else
    disp('on rejette l hypothèse que coefnwifeinc = 0.01 à 5%');
end

disp('p-value : ')
disp(tdis_prb(estim10,n-k))


%% Question 11

SSR110 = u8'*u8;

[n,k]=size(mroz);
X11 = [ones(n,1) educ exper kidslt6 kidsge6];
y11 = lwage-0.01.*nwifeinc-0.05.*city;
beta11 = inv(X11'*X11)*X11'*y11;
u11 = y11-X11*beta11;
SSR111 = u11'*u11;

[n,k]=size(X8);
F11 = ((SSR111-SSR110)/SSR110)*(n-k)/2;

if fdis_prb(F11,2,n-k)>=0.05
    disp('on accepte l hypothèse que coefnwifeinc = 0.01 et que coefcity = 0.05 à 5%');
else
    disp('on rejette l hypothèse que coefnwifeinc = 0.01 et que coefcity = 0.05 à 5%');
end

disp('p-value : ')
disp(fdis_prb(F11,2,n-k))


%% Question 12

import gridfit.*;

X = educ;
Y = exper;
Z = wage;
idx = any(isnan([X Y Z]),2);
X(idx) = []; 
Y(idx) = []; 
Z(idx) = [];
[~,idx] = unique([X Y Z], 'rows');
X = X(idx); 
Y = Y(idx); 
Z = Z(idx);
xi = linspace(min(X), max(X), 50);
yi = linspace(min(Y), max(Y), 50);
zi = gridfit(X, Y, Z, xi, yi);
figure('Name','Modélisation 2D de wage en fonction de exper et educ','NumberTitle','off');
plot3(X,Y,Z,'r*')
xlabel('educ')
ylabel('exper')
hold on
surf(xi, yi, zi)
axis tight vis3d
grid on


%% Question 13

kids = kidsge6-kidslt6;
X13 = [ones(n,1) city educ exper nwifeinc kids];
y13 = wage;
beta13 = inv(X13'*X13)*X13'*y13;
u13 = y13-X13*beta13;

[n,k] = size(X13);
sig13 = u13'*u13/(n-k);
std13 = sqrt(diag(sig13*inv(X13'*X13)));
t13 = beta13./std13;
estim13 = t13(6);

if estim9<=tdis_inv(0.975,n-k) && estim9>=tdis_inv(0.025,n-k)
    disp('on accepte l hypothèse que kidsge6 = kidslt6 à 5%');
else
    disp('on rejette l hypothèse que kidsge6 = kidslt6 à 5%');
end

disp('p-value : ')
disp(tdis_prb(estim13,n-k))

%% Question 14

y14 = u7.^2;
X14 = X7;
[n,k] = size(X14);
beta14 = inv(X14'*X14)*X14'*y14;
u14 = y14-X14*beta14;
sig14 = u14'*u14/(n-k);
std14 = sqrt(diag(sig14*inv(X14'*X14)));
t14 = (beta14)./std14;
SSR140 = u14'*u14;

X141 = X14(:,[1]);
y141 = y14;
beta141 = inv(X141'*X141)*X141'*y141;
u141 = y141-X141*beta141;
SSR141 = u141'*u141;

F14 = ((SSR141-SSR140)/SSR140)*(n-k)/6;
disp('p-value : ');
disp(fdis_prb(F14,6,n-k));



val14cor = sqrt(kidslt6)+1;
[n,k] = size(mroz);
X14cor = X14./val14cor;
y14cor = y14./val14cor;
[n,k] = size(X14cor);
beta14cor = inv(X14cor'*X14cor)*X14cor'*y14cor;
u14cor = y14cor-X14cor*beta14cor;

y14cor = u14cor.^2;
beta14cor = inv(X14cor'*X14cor)*X14cor'*y14cor;
u14cor = y14cor-X14cor*beta14cor;
SSR14cor0 = u14cor'*u14cor;

X14cor1 = X14cor(:,[1]);
y14cor1 = y14cor;
beta14cor1 = inv(X14cor1'*X14cor1)*X14cor1'*y14cor1;
u14cor1 = y14cor1-X14cor1*beta14cor1;
SSR14cor1 = u14cor1'*u14cor1;

F14cor = ((SSR14cor1-SSR14cor0)/SSR14cor0)*(n-k)/6;
disp('p-valuecor : ');
disp(fdis_prb(F14cor,6,n-k));



%% Question 15

womenu43 = mroz(:,5)<=43;
womena43 = 1-womenu43;
[n,k]=size(mroz);

X15u43 = [ones(n,1) educ.*womenu43 exper.*womenu43 city.*womenu43 nwifeinc.*womenu43 kidslt6.*womenu43 kidsge6.*womenu43];
y15u43 = lwage.*womenu43;
beta15u43 = inv(X15u43'*X15u43)*X15u43'*y15u43;
u15u43 = y15u43-X15u43*beta15u43;
SSR15u43 = u15u43'*u15u43;

X15a43 = [ones(n,1) educ.*womena43 city.*womena43 exper.*womena43 nwifeinc.*womena43 kidslt6.*womena43 kidsge6.*womena43];
y15a43 = lwage.*womena43;
beta15a43 = inv(X15a43'*X15a43)*X15a43'*y15a43;
u15a43 = y15a43-X15a43*beta15a43;
SSR15a43 = u15a43'*u15a43;


X15 = X8;
y15 = y8;
[n,k] = size(X15);
beta15 = inv(X15'*X15)*X15'*y15;
u15 = y15-X15*beta15;
SSR15 = u15'*u15;

F15=(SSR15-(SSR15a43+SSR15u43))*(n-k)/(7*(SSR15a43+SSR15u43));
disp('p-value : ');
disp(fdis_prb(F15,6,n-k));



womenu30 = mroz(:,5)<=30;
womenb3043 = 1-(womena43+womenu30);
[n,k]=size(mroz);
X15bis = [ones(n,1) educ.*womenu30 exper.*womenu30 city.*womenu30 nwifeinc.*womenu30 kidslt6.*womenu30 kidsge6.*womenu30 educ.*womenb3043 exper.*womenb3043 city.*womenb3043 nwifeinc.*womenb3043 kidslt6.*womenb3043 kidsge6.*womenb3043 educ.*womena43 city.*womena43 exper.*womena43 nwifeinc.*womena43 kidslt6.*womena43 kidsge6.*womena43];
y15bis = y15;
beta15bis = inv(X15bis'*X15bis)*X15bis'*y15bis;
u15bis= y15bis-X15bis*beta15bis;
[n,k]=size(X15bis);
sig15bis = u15bis'*u15bis/(n-k);
std15bis = sqrt(diag(sig15bis*inv(X15bis'*X15bis)));
t15bis =(beta15bis)./std15bis;
disp('p-values : ');
disp(tdis_prb(t15bis,n-k));



%% Question 16

X16 = [womenu30 womenb3043 womena43];
y16 = y8;

beta16 = inv(X16'*X16)*X16'*y16;
u16 = y16-X16*beta16;
sig16 = u16'*u16/(n-k);
std16 = sqrt(diag(sig16*inv(X16'*X16)));
t16 = (beta16)./std16;

disp('p-values : ');
disp(tdis_prb(t16,n-k));