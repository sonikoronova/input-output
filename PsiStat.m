function res=PsiStat(x, xt)
    xt(xt==0)=1e-2;
    x(x==0)=1e-2;
    s=0.5*(abs(x)+abs(xt));
    res=sum(sum(abs(xt).*abs(log(abs(xt./s)))+abs(x).*abs(log(abs(x./s)))));
    if isnan(res)
        display("SHEET")
    end
    res=res/sum(sum(xt));
end