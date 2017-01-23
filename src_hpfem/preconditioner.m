function [out]=preconditioner(in,arg)
Xzz=arg.Xzz;
nz=arg.nz;

XEgEg=arg.XEgEg;
nheg=arg.nheg;

XFgFg=arg.XFgFg;
nhfg=arg.nhfg;

XFF=arg.XFF;
nhf=arg.nhf;

XIgIg=arg.XIgIg;
nhig=arg.nhig;

XII=arg.XII;
nhi=arg.nhi;

% Zero block
inz=in(1:nz);

tmp=arg.L\(arg.P*inz);
tmp2=arg.U\tmp;
outz=arg.Q*tmp2;
%outz=Xzz\inz;

%Edge block
if nheg > 0
    inEg=in(nz+1:nz+nheg);
    outEg=pcg(XEgEg,inEg,1e-8,400);
else
    outEg=[];
end
%Face block
if nhfg > 0
    inFg=in(nz+nheg+1:nz+nheg+nhfg);
    outFg=pcg(XFgFg,inFg,1e-8,400);
else
    outFg=[];
end

%Face block (non-gradients)
if nhf > 0
    inF=in(nz+nheg+nhfg+1:nz+nheg+nhfg+nhf);
    %outF=XFF\inF;
    tmpF=arg.Lf\(arg.Pf*inF);
    tmpF2=arg.Uf\tmpF;
    outF=arg.Qf*tmpF2;
else
    outF=[];
end

%Interior block (gradients)
if nhig > 0
    inIg=in(nz+nheg+nhfg+nhf+1:nz+nheg+nhfg+nhf+nhig);
    outIg=pcg(XIgIg,inIg,1e-8,400);
else
    outIg=[];
end

%Interior block (non-gradients)
if nhi > 0
    inI=in(nz+nheg+nhfg+nhf+nhig+1:nz+nheg+nhfg+nhf+nhig+nhi);
    %outI=XII\inI;
    tmpI=arg.Li\(arg.Pi*inI);
    tmpI2=arg.Ui\tmpI;
    outI=arg.Qi*tmpI2;
else
    outI=[];
end


out=[outz;outEg;outFg;outF;outIg;outI];
