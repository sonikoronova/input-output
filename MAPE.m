function res=MAPE(x, xt)
    [n,m]=size(x);
    xt(xt==0)=1e-2;
    res=100/(n*m)*sum(sum(abs(x-xt)./abs(xt)));
end