% COMPUTATIONAL EFFORT OF MMSE
NT=3:8;
NR=NT;

% MMSE-SQRD
f_SQRD=4/3*NT.^3+4*NT.^2.*NR+1/3*NT.^2+2*NT.*NR+25/6*NT;
plot(NT,f_SQRD/1000,'o-k');
hold on;

% MMSE-QR 
f_QR=f_SQRD-2*NT.^2-2*NT;
plot(NT,f_QR/1000,'*-k');
hold on;

%MMSE-SQRD-PSA
f_PSA=14/3*NT.^3+4*NT.^2.*NR+27/2*NT.^2+3*NT.*NR+89/6*NT-7*NR-30;
plot(NT,f_PSA/1000,'d-k');
hold on;

%WORST CASE
f_WORST=f_SQRD+f_PSA;
plot(NT,f_WORST/1000,'k');
hold on;

xlabel('NT=NR');
ylabel('10^3 flops');
title('Number of operation f in flops for unsorted MMSE-QRD,MMSE-SQRD,MMSE-SQRD-PSA and the algorithm by Hassibi.')
legend('MMSE-SQRD','MMSE-QR','MMSE-SQRD-PSA','worst case');