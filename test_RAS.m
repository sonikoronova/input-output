% данные взяты из статьи "A short note on RAS method" 
% http://www.scienpress.com/Upload/AMAE/Vol%203_4_12.pdf

format short;

A0=[21.2, 29.0, 19.9; 
    9.0, 48.8, 22.4;
    49.8, 62.2, 7.8];

u=[75; 82; 125];

v=[85, 150, 45];

% вызовем наш RAS
%A=GRAS(A0, u, v);
%A=MyGRAS_1(A0, u, v);
A=INSD(A0, u, v);

%err
%A
%sum(A, 1)
%sum(A, 2)

MAPE(A, A0)
WAPE(A, A0)
SWAD(A, A0)
PsiStat(A, A0)
RSQ(A, A0)