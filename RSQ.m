function res=RSQ(x, xt)
    xt(xt==0)=1e-6;
    res=corr2(x,xt)*corr2(x,xt);
end