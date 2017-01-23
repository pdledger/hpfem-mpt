function [out]=preconditioner_dir(in,arg)
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

%Edge block
if nheg > 0
    inEg=in(nz+1:nz+nheg);
    tmpE=arg.Leg\(arg.Peg*inEg);
    tmpE2=arg.Ueg\tmpE;
    outEg=arg.Qeg*tmpE2;
else
    outEg=[];
end
%Face block
if nhfg > 0
    inFg=in(nz+nheg+1:nz+nheg+nhfg);
    tmpF=arg.Lfg\(arg.Pfg*inFg);
    tmpF2=arg.Ufg\tmpF;
    outFg=arg.Qfg*tmpF2;
else
    outFg=[];
end

%Face block (non-gradients)
if nhf > 0
    inF=in(nz+nheg+nhfg+1:nz+nheg+nhfg+nhf);
    tmpF=arg.Lf\(arg.Pf*inF);
    tmpF2=arg.Uf\tmpF;
    outF=arg.Qf*tmpF2;
else
    outF=[];
end

%Interior block (gradients)
if nhig > 0
    inIg=in(nz+nheg+nhfg+nhf+1:nz+nheg+nhfg+nhf+nhig);
    tmpI=arg.Lig\(arg.Pig*inIg);
    tmpI2=arg.Uig\tmpI;
    outIg=arg.Qig*tmpI2;
else
    outIg=[];
end

%Interior block (non-gradients)
if nhi > 0
    inI=in(nz+nheg+nhfg+nhf+nhig+1:nz+nheg+nhfg+nhf+nhig+nhi);
    tmpI=arg.Li\(arg.Pi*inI);
    tmpI2=arg.Ui\tmpI;
    outI=arg.Qi*tmpI2;
else
    outI=[];
end

out=[outz;outEg;outFg;outF;outIg;outI];