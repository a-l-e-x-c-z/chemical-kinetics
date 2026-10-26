% kinetikaBZ

% volby pro nastaveni
provest_pokusny_vypocet = false; % true, false
provest_fitovani = true;
detailni_vypis_fitovani = true;
vypis_fitovani_hvezdicky = false; % misto detailniho vypisuje postupne hvezdicky

% prvni nastrel pro rychlostni koeficienty
k1 = 0.05       % kAB, A ---> B
k2 = 0.001          % kBA, B ---> A
k3 = 0.005       % kBC, B ---> C
k4 = 0.005       % kCD, C ---> D
k5 = 0.005          % kBE, B ---> E
k6 = 0.005          % kEF, E ---> F
k7 = 0.005          % kAZ, A ---> Z

% vektor rychlostnich koeficientu k1, k2, k3, k4........
koeficienty = [k1 k2 k3 k4 k5 k6 k7];

% rozsah/omezeni hodnot kinetickych koeficientu v seznamu
% pokud se ve vypoctu rozsah prekroci, tak se koeficient prenastavi
% aby byl opet v rozsahu
koeficienty_rozsah = [0.0  100.0];

% funkce pro soustavu kinetickych rovnic
%  t = casova osa
%  koncentrace = vektor koncentraci slozek A, B, C, D.....
%   koeficienty = vektor kinetickych koeficienty kAB, kAC, kBC.....
%
function dkonc_dt = diff_rovnice(t, koncentrace, koeficienty)

% koncentrace jednotlivych slozek A, B, C, D....
% jsou ulozeny ve vektoru koncentrace
A = koncentrace(1); % prvni polozka vektoru
B = koncentrace(2); % druha polozka vektoru
C = koncentrace(3);
D = koncentrace(4);
E = koncentrace(5);
F = koncentrace(6);
Z = koncentrace(7);
%D = koncentrace(4);

% kineticke koeficienty kAB, kBC, kCB, kCD....
% jsou ulozeny ve vektoru koeficienty
kAB = koeficienty(1);
kBA = koeficienty(2);
kBC = koeficienty(3);
kCD = koeficienty(4);
kBE = koeficienty(5);
kEF = koeficienty(6);
kAZ = koeficienty(7);
%kCD = koeficienty(5);

% zapis chemickych/kinetickych rovnic s diferencemi, koncentracemi a koeficienty
% A ---kAB---> B     kAB = k1
% B ---kBA---> A     kBA = k2
% B ---kBC---> C     kBC = k3
% C ---kCD---> D     kCD = k4
% B ---kBE---> E     kBE = k5
% E ---kEF---> F     kEF = k6
% A ---kAZ---> Z     kEF = k7
dAdt = -kAB*A -kAZ*A+ kBA*B;
dBdt = -kBA*B - kBC*B - kBE*B + kAB*A;
dCdt = -kCD*C + kBC*B;
dDdt = kCD*C;
dEdt = -kEF*E + kBE*B;
dFdt = kEF*E;
dZdt = kAZ*A;
% dDdt = kCD*C;

% zapsani vypoctenych diferenci do vystupniho vektoru
dkonc_dt(1) = dAdt;
dkonc_dt(2) = dBdt;
dkonc_dt(3) = dCdt;
dkonc_dt(4) = dDdt;
dkonc_dt(5) = dEdt;
dkonc_dt(6) = dFdt;
dkonc_dt(7) = dZdt;
%dkonc_dt(4) = dDdt;

endfunction

% rozsah casu pro vypocet diferencialnich rovnic, t=0 az t=tmax
tmax = 150.0;

% experimentalni data casy a koncentrace A, B, C... ve vektoru
% prvni sloupecek jsou rostouci casy t
% dalsi sloupecek je koncentrace A, dalsi sloupecky koncentraci
% cas t    koncentrace A    koncentrace B    koncentrace C
t_koncentrace_experiment = [
0.0     1.0   0.0   0.0   0.0   0.0   0.0 0.0;   % je vhodne, aby prvni radek odpovidal casu t=0
28.1294	0.46103	0.19045	0.04924	0.00854	0.06598	0.00282	0.22195;
28.1294	0.52064	0.16103	0.06571	0.01012	0.08537	0.00438	0.15276;
42.19409	0.27648	0.19846	0.08699	0.02221	0.12453	0.00866	0.28267;
42.19409	0.34242	0.23293	0.08208	0.01843	0.12647	0.00894	0.18873;
56.25879	0.15864	0.16456	0.11694	0.04048	0.18896	0.01748	0.31294;
56.25879	0.20448	0.19822	0.11785	0.03509	0.19464	0.01837	0.23135;
]

% popisky kinetickych grafu
titulek = "Kineticke rovnice";
popis_x = "osa x, cas";
popis_y = "osa y, koncentrace A, B, C, D, E, F";

% nasledujici funkce provadeji vlastni vypocty a neni je vetsinou treba menit

% objektivni funkce pro hledani minimalni diference vypoctenych
% a experimentalnich dat
% pomoci resitele diferencialnich rovnic ode23 se pocitaji hodnoty koncentraci
function of = objektivni_funkce(koeficienty, koeficienty_rozsah, tkonc_exp, log, stars);

soucet_diff = 0;

% koeficienty musi byt v rozsahu
for i=1:length(koeficienty)
  if koeficienty(i)<koeficienty_rozsah(1) || koeficienty(i)>koeficienty_rozsah(2)
    koeficienty(i) = (koeficienty_rozsah(1) + koeficienty_rozsah(2))/2.0;
  endif
end

% iterace prez vsechny experimentalni casy
pocet_experimentu = size(tkonc_exp)(1);
for i=2:pocet_experimentu
   t_span = [0.0 tkonc_exp(i,1)];
   konc_start = tkonc_exp(1,2:size(tkonc_exp)(2));

   % kineticke rovnice ve tvaru pro resitele ODE
   ode_rovnice = @(t,koncentrace)diff_rovnice(t,koncentrace, koeficienty);

   % resitel sady diferencialnich rovnic
   [t_vyp, konc_vyp] = ode23(ode_rovnice, t_span, konc_start);

   tvpm = length(t_vyp);

   % rozdil koncentraci v experimentalnich a vypoctenych hodnotach
   diff_A = tkonc_exp(i,2)-konc_vyp(tvpm,1);
   diff_B = tkonc_exp(i,3)-konc_vyp(tvpm,2);
   diff_C = tkonc_exp(i,4)-konc_vyp(tvpm,3);
   diff_D = tkonc_exp(i,5)-konc_vyp(tvpm,4);
   diff_E = tkonc_exp(i,6)-konc_vyp(tvpm,5);
   diff_F = tkonc_exp(i,7)-konc_vyp(tvpm,6);
   diff_Z = tkonc_exp(i,8)-konc_vyp(tvpm,7);
   %diff_D = tkonc_exp(i,5)-konc_vyp(tvpm,4);

   % objektivni funkce jako soucet rozdilu vypoctenych a experimentalnich koncentraci
   soucet_diff = soucet_diff + (diff_A*diff_A + diff_B*diff_B + diff_C*diff_C + diff_D*diff_D + diff_E*diff_E + diff_F*diff_F + diff_Z*diff_Z)/7.0;

   if (log == true)
       disp("iterace t,cas: "), disp(tkonc_exp(i,1)), disp("koncentrace: "), disp(tkonc_exp(i,2:end)), disp(konc_vyp(tvpm));
   endif
endfor
of = soucet_diff;

if (log == true)
  disp("of= "), disp(soucet_diff), disp(koeficienty);
endif

if (stars == true)
  fprintf("*");
endif

endfunction

% objektivni funkce ve tvaru pro fminsearch
% fminsearch je metoda pro hledani minima funkce
of_rovnice = @(koeficienty)objektivni_funkce(koeficienty, koeficienty_rozsah, t_koncentrace_experiment, detailni_vypis_fitovani, vypis_fitovani_hvezdicky);

% vlastni hledani minima objektivni funkce na zaklade kinetickych koeficientu
if (provest_fitovani == true)
 vystupni_koeficienty = fminsearch(of_rovnice, koeficienty);

 fprintf("vystupni koeficienty k1, k2, k3.....\n");
 vystupni_koeficienty
 fprintf("-----------------------------------------------\n");

% vypsani vysledku a vykresleni grafu kinetiky
 casovy_rozsah = [0 tmax];
 %vystupni_koeficienty = [0.93 0.89 0.87 0.88];

 vypocet = @(t,koncentrace)diff_rovnice(t,koncentrace, vystupni_koeficienty);
 [t_vypocet, konc_vypocet] = ode23(vypocet, casovy_rozsah, t_koncentrace_experiment(1,2:size(t_koncentrace_experiment)(2)));

 plot(t_vypocet,konc_vypocet,t_koncentrace_experiment(:,1), t_koncentrace_experiment(:,2:size(t_koncentrace_experiment)(2)),'o');
 xlim([0 tmax]);
 title(titulek);
 xlabel(popis_x);
 ylabel(popis_y);

 %vypis vypoctenych dat
 fprintf("vypoctene kinetiky: cas koncentrace...... \n");
 for i=1:size(t_vypocet)(1)
  fprintf("%f ", t_vypocet(i));
  for j=1:size(konc_vypocet)(2)
    fprintf("%f ", konc_vypocet(i,j));
  endfor
  fprintf("\n");
endfor

endif

% pomocny vypocet kinetiky
% pomocna sada koeficientu a hodnot
if provest_pokusny_vypocet == true
pokusne_koeficienty = [0.05 0.01 0.05 0.05 0.05 0.05 0.05]; % kineticke koeficienty k1, k2, k3.....
pokusny_casovy_rozsah = [0.0 150.0];    % cas 0 az cas 60
pokusne_pocatecni_koncentrace = [1.0 0.0 0.0 0.0 0.0 0.0 0.0]; % cas 0, koncentrace A, B, C

% diff_rovnice(t, koncentrace, koeficienty)
pokusny_vypocet = @(t,koncentrace)diff_rovnice(t,koncentrace, pokusne_koeficienty);
[t_pokusne, konc_pokusne] = ode23(pokusny_vypocet, pokusny_casovy_rozsah, pokusne_pocatecni_koncentrace);

fprintf("pokusny vypocet kinetiky: cas koncentrace...... \n");
for i=1:size(t_pokusne)(1)
  fprintf("%f ", t_pokusne(i));
  for j=1:size(konc_pokusne)(2)
    fprintf("%f ", konc_pokusne(i,j));
  endfor
  fprintf("\n");
endfor


plot(t_pokusne,konc_pokusne,t_pokusne, konc_pokusne,'o');
xlim(pokusny_casovy_rozsah);
title(titulek);
xlabel(popis_x);
ylabel(popis_y);
endif

