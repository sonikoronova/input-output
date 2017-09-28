function res=SWAD(x, xt)
    xt(xt==0)=1e-6;
    res=sum(sum(abs(xt).*abs(x-xt)))/sum(sum(xt.*xt));
end