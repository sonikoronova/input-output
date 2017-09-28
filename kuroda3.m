function A=kuroda3(A0, u,v)

    [n,m]=size(A0);

    A0(A0==0)=1e-6;
    u0=sum(A0, 2);
    v0=sum(A0, 1);
    coefs=zeros(1, n*m);
    c=zeros(1, n*m);
    for i=1:n
        for j=1:m
            coefs(m*(i-1)+j)=1/(u(i)*u(i))+1/(v(j)*v(j));
            c(m*(i-1)+j)=-A0(i,j)*( 1/(u(i)*u0(i))+1/(v(j)*v0(j)) );
        end
    end
    Q=0.5*diag(coefs);

    A=zeros(n+m, n*m);
    for i=1:n
        for j=(m*(i-1)+1):(m*i)
            A(i,j)=1;
        end
    end
    for i=1:m
        for j=i:n:(n*m)
            A(n+i,j)=1;
        end
    end
    rhs=[u' v];

    clear model;
    model.Q=sparse(Q);
    model.A=sparse(A);
    model.rhs=rhs;
    model.obj=c;
    model.sense='<';

    results = gurobi(model);
    A=reshape(results.x, [n,m])';
end