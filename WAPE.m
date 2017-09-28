function res=WAPE(x, xt)
    xt(xt==0)=1e-6;
    res=100/sum(sum(abs(xt)))*sum(sum(abs(x-xt)));
end