function hcurlerr =hcurl(nelem,gorder,eltype,edgecof,facecof,order,esizet,...
    intma,glob,globfa,unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,sol,coord,...
    nip,intxi,inteta,intzeta,intw,ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,...
    ecphy1,ecphz1,ecphx2,ecphy2,ecphz2,probdata,gphx1,gphy1,gphz1,...
    gphx2,gphy2,gphz2,phh11,phh12)


lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

error=0;
errord=0;


for i=1:nelem
    
    xy = coord(intma(i,1:4),1:3);
    
    flag=0;
    for j=1:6
        for pp=1:gorder
            for k=1:3
                lec(j,pp,k)=edgecof(glob(i,j),((pp-1)*3)+k);
                if(abs(edgecof(glob(i,j),((pp-1)*3)+k))>0.00000001)
                    flag=1;
                end
            end
        end
    end
    
    for j=1:4
        for pp=1:gorder*(gorder-1)/2
            for k=1:3
                lfc(j,pp,k)=facecof(globfa(i,j),((pp-1)*3)+k);
            end
        end
    end
    
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
    
    mycoord = zeros(gesizet,3);
    
    %transfer coefficents to locations vertices
    mycoord(1:4,1:3) = xy(1:4,1:3);
    
    % edge functions
    for ii=1:6
        for p=1:gorder
            for j=1:3
                mycoord(4+ii+6*(p-1),j)=lec(ii,p,j);
            end
        end
    end
    
    % face functions
    for iii=1:4
        for ii=1:(gorder-1)*gorder/2
            for j=1:3
                mycoord(4+6*gorder+(iii-1)*gorder*(gorder-1)/2+ii,j)= lfc(iii,ii,j);
            end
        end
    end
    
    if flag==0
        gesizet=4;
    end
    
    %-------------------------------------------------------------------------
    % work out numbering of basis functions
    bhelp=rowfun(unkz,unksid,unkint,order,glob,globfa,...
        i,unkfatp1,unkfatp2,unkfatp3,esizet);
    
    %         % extract local solution for this element
    lsol=zeros(esizet,1);
    for ii=1:esizet
        row=bhelp(ii);
        if row>0
            lsol(ii,1)=sol(row,1);
        end
    end
    
    % choose correct set of stored basis functions
    if eltype(i)==1
        gphx=gphx1;
        gphy=gphy1;
        gphz=gphz1;
        phh1=phh11;
        ephx=ephx1;
        ephy=ephy1;
        ephz=ephz1;
        ecphx=ecphx1;
        ecphy=ecphy1;
        ecphz=ecphz1;
    else
        gphx=gphx2;
        gphy=gphy2;
        gphz=gphz2;
        phh1=phh12;
        ephx=ephx2;
        ephy=ephy2;
        ephz=ephz2;
        ecphx=ecphx2;
        ecphy=ecphy2;
        ecphz=ecphz2;
    end
    
    
    %-------------------------------------------------------------------------
    % evaluate covairant mapping (for linear geometry mapping is constant)
    if flag==0
        gph(1:gesizet,1)=gphx(1,1:gesizet)';
        gph(1:gesizet,2)=gphy(1,1:gesizet)';
        gph(1:gesizet,3)=gphz(1,1:gesizet)';
        
        [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        
        for pp=1:nip
            
            ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
            % use stored functions
            ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                (ephy(pp,1:esizet)'*aeta(1:3))+...
                (ephz(pp,1:esizet)'*azeta(1:3)); 
            
            
            % compute the solution for this problem, at this integration point
            % in this element
            e=zeros(3,1);
            e=ph'*lsol;
            
            % use stored functions
            ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                (ecphy(pp,1:esizet)'*aseta(1:3))+...
                (ecphz(pp,1:esizet)'*aszeta(1:3));
            
            % compute the solution for this problem, at this integration point
            % in this element
            curle=zeros(3,1);
            curle=ph'*lsol;
            
            
            % computx x,y,z
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            
            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg);                        
            
            fun=probdata.es.exactcurlfun;
            arg=probdata.es.exactcurlfunarg;
            [curlexact,dum]=fun(x,y,z,arg);
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(e(ii)-exact(ii))^2+...
                    abs(curle(ii)-curlexact(ii))^2);
                errord=errord+intw(pp)*det*(abs(exact(ii))^2+abs(curlexact(ii))^2);
            end
            
            
        end
        
    else
        
        for pp=1:nip
            gph(1:gesizet,1)=gphx(pp,1:gesizet)';
            gph(1:gesizet,2)=gphy(pp,1:gesizet)';
            gph(1:gesizet,3)=gphz(pp,1:gesizet)';
            ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
            
            % evaluate covairant mapping
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
            
            %             % use stored functions            
            ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                (ephy(pp,1:esizet)'*aeta(1:3))+...
                (ephz(pp,1:esizet)'*azeta(1:3));
                  
            % compute the solution for this problem, at this integration point
            % in this element
            e=zeros(3,1);            
            e=ph'*lsol;
            
            
            %             % use stored functions
            ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                (ecphy(pp,1:esizet)'*aseta(1:3))+...
                (ecphz(pp,1:esizet)'*aszeta(1:3));
            
            % compute the solution for this problem, at this integration point
            % in this element
            curle=zeros(3,1);
            curle=ph'*lsol;
            
            
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg);
            
            
            
            fun=probdata.es.exactcurlfun;
            arg=probdata.es.exactcurlfunarg;
            [curlexact,dum]=fun(x,y,z,arg);
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(e(ii)-exact(ii))^2+...
                    abs(curle(ii)-curlexact(ii))^2);
                errord=errord+intw(pp)*det*(abs(exact(ii))^2+abs(curlexact(ii))^2);
            end
            
        end
        
    end
    
    % end of loop over elements
end
hcurlerr = sqrt(error/errord);
display(['The hcurl error is = ',num2str(hcurlerr)])