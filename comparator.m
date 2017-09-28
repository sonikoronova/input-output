A0=xlsread('ukrstat2015.xlsx', 'Лист1', 'B6:T24');
At=xlsread('ukrstat2014.xlsx', 'Sheet1', 'B6:T24');
u=xlsread('ukrstat2014.xlsx', 'Sheet1', 'U6:U24');
v=xlsread('ukrstat2014.xlsx', 'Sheet1', 'B25:T25');

A0(isnan(A0))=1e-2;
At(isnan(At))=1e-2;
A=GRAS(A0, u, v);

mape=MAPE(A, At)
wape=WAPE(A, At)
swad=SWAD(A, At)
psi=PsiStat(A, At)
rsq=RSQ(A, At)