function [Afin,r,s] = MyGRAS_1(A,u,v,sq,epsl)
% PURPOSE: estimate a new matrix Afin on the base of A, u and v. Also, sq should
% be 2, if there are a lot of zeros in A, otherwise sq = 1 automaticaly.
% -------------------------------------------------------------------------
% USAGE: Afin = MyGRAS(A,u,v) OR [Afin,r,s] = MyGRAS(A,u,v,2) with or
% without epsl included as the fourth argument, where
% INPUT:
% -> A = benchmark (base) matrix, not necessarily square
% -> u = needed row sums, size:  n x 1
% -> v = needed column sums, size:  1 x m
% -> epsl = convergence tolerance level; if empty, the default threshold
% is 0.1e-5 (=0.000001)
% OUTPUT:
% -> Afin = estimated/adjusted/updated matrix
% -> r = substitution effects (row multipliers)
% -> s = fabrication effects (column multipliers)
% -------------------------------------------------------------------------
% REFERENCES: 1) Junius T. and J. Oosterhaven (2003), The solution of
% updating or regionalizing a matrix with both positive and negative
% entries, Economic Systems Research, 15, pp. 87-96.
% 2) Lenzen M., R. Wood and B. Gallego (2007), Some comments on the GRAS
% method, Economic Systems Research, 19, pp. 461-465.
% 3) Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note on the
% GRAS method, Economic Systems Research, XX, pp. XX-XX.
% -------------------------------------------------------------------------
% Written by: Denis Sokolov, 9/05/2015
% Current e-mail: sokol0077@gmail.com

[n,m] = size (A);

%for avoiding zeros we replace them with 0.1e-8, in the end we will replace
%them back
for i=1:1:n
    if(u(i)==0)
        u(i)=0.1e-8;
    end;
end;

for j=1:1:m
    if(v(j)==0)
        v(j)=0.1e-8;
    end;
end;

format long;

%creating an initial vector of coefficients r (for rows)
for i=1:1:n
    r(i,1)=1;
end;

%creating an initial vector of coefficients s (for columns)
for i=1:1:m
    s(1,i)=1;
end;

Iter_number = 0; %total number of iterrations

if nargin < 5 || isempty(epsl)
    epsl = 0.1e-5; %default tolerance level
end;

if nargin < 5 || isempty(sq)
    sq = 1; %default sq
end;

otkl = 10;

%GRAS

%Constructing P and N
P=zeros(n,m);
P(A>0)=A(A>0);
N=P-A;

if(sq==1)
    while (otkl>epsl)%while deviation is more than eps.
        Iter_number = Iter_number+1;%one new iterration

        for j=1:1:n
            summStrP=0;
            summN=0;
            for k=1:1:m
                summStrP=summStrP+(P(j,k)*s(1,k));
                summN=summN+(N(j,k)/s(1,k));
            end;
            if (summStrP>0)
                r(j,1)=((u(j,1)+sqrt((u(j,1).^2 + 4*summStrP*summN)))...
                    /(2*summStrP));
            else
                r(j,1)=-1*(summN/u(j,1));
            end;
        end;
        for k=1:1:m
            summStlbP=0;
            summN=0;
            for j=1:1:n
                summStlbP=summStlbP+(P(j,k)*r(j,1));
                summN=summN+(N(j,k)/r(j,1));
            end;
            if (summStlbP>0)
                s(1,k)=((v(1,k)+sqrt((v(1,k).^2 + 4*summStlbP*summN)))...
                    /(2*summStlbP));
            else
                s(1,k)=-1*(summN/v(1,k));
            end;
        end;

        %calculating Aprom matrix (that we have right now) for comparing with
        %u and v
        for i=1:1:n
            for j=1:1:m
                if A(i,j)<0
                    Aprom(i,j)=(1./r(i,1))*A(i,j)*(1./s(1,j));
                else
                    Aprom(i,j)=(r(i,1)*A(i,j)*s(1,j));
                end;
            end;
        end;

        %calculating transitional flanking rows sum
        for i=1:1:n
            uprom(i,1)=0;
            for j=1:1:m
                uprom(i,1)=uprom(i,1)+Aprom(i,j);
            end;
        end;

        %calculating transitional flanking columns sum
        for j=1:1:m
            vprom(1,j)=0;
            for i=1:1:n
                vprom(1,j)=vprom(1,j)+Aprom(i,j);
            end;
        end;

        %calculating deviation in flanking rows sum
        for i=1:1:n
            urasn(i,1)=abs(u(i,1)-uprom(i,1));
        end;

        %calculating deviation in flanking columns sum
        for j=1:1:m
            vrasn(1,j)=abs(v(1,j)-vprom(1,j));
        end;

        otkl=max(max(urasn(:,1)),max(vrasn(1,:)));%the biggest отклонение
    end;
else
    for g = 1:1:100
        Iter_number = Iter_number+1;%one new iterration

        for j=1:1:n
            summStrP=0;
            summN=0;
            for k=1:1:m
                summStrP=summStrP+(P(j,k)*s(1,k));
                summN=summN+(N(j,k)/s(1,k));
            end;
            if (summStrP>0)
                r(j,1)=((u(j,1)+sqrt((u(j,1).^2 + 4*summStrP*summN)))...
                    /(2*summStrP));
            else
                r(j,1)=-1*(summN/u(j,1));
            end;
        end;
        for k=1:1:m
            summStlbP=0;
            summN=0;
            for j=1:1:n
                summStlbP=summStlbP+(P(j,k)*r(j,1));
                summN=summN+(N(j,k)/r(j,1));
            end;
            if (summStlbP>0)
                s(1,k)=((v(1,k)+sqrt((v(1,k).^2 + 4*summStlbP*summN)))...
                    /(2*summStlbP));
            else
                s(1,k)=-1*(summN/v(1,k));
            end;
        end;
    end;
end;


%creating the final matrix Afin
for i=1:1:n
    for j=1:1:m
        if A(i,j)<0
            Afin(i,j)=(1./r(i,1))*A(i,j)*(1./s(1,j));
        end;
        if A(i,j)>0.1e-6
            Afin(i,j)=(r(i,1)*A(i,j)*s(1,j));
        end;
        if A(i,j)>=0 && A(i,j)<0.1e-6
            Afin(i,j)=0;
        end;
    end;
end;










