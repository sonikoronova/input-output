function A=INSD(A0, u, v)
    [n,m]=size(A0);

    A0(A0==0)=1e-6;

    coefs=reshape(A0', [], n*m);
    Q=diag(1./coefs);
    c=-2*ones(1, n*m);

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