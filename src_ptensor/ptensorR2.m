% Function to compute Rank 2 polarization tensor

function ptensor =ptensorR2(nelem,gorder,eltype,edgecof,facecof,order,esizet,...
    intma,glob,globfa,sol,coord,unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,...
    nip,intxi,inteta,intzeta,intw,ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,...
    ecphy1,ecphz1,ecphx2,ecphy2,ecphz2,probdata,matc,mat,gphx1,gphy1,gphz1,...
    gphx2,gphy2,gphz2,phh11,phh12);


lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

delta=probdata.matr.delta;
mu=probdata.matr.mu;
muz=probdata.matr.muz;
omega=probdata.matr.omega;
sigma=probdata.matr.sigma;
nrhs=probdata.jb.nrhs;

cchkrank2=zeros(3,3);
nrank2=zeros(3,3);

term1=0;
term2=zeros(3,3);


myT=[];
for i=1:nelem
    if matc(i)==1
        
        xy = coord(intma(i,1:4),1:3);
        myT=[myT; intma(i,:)];
        
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
        lsol=zeros(esizet,3);
        for ii=1:esizet
            row=bhelp(ii);
            if row>0
                lsol(ii,1:3)=conj(sol(row,1:3));
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
                e=zeros(3,3);
                e=ph'*lsol;
                                
                
                % use stored functions
                ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                        (ecphy(pp,1:esizet)'*aseta(1:3))+...
                        (ecphz(pp,1:esizet)'*aszeta(1:3));
                
                % compute the solution for this problem, at this integration point
                % in this element
                curle=zeros(3,3);
                curle=ph'*lsol;
                
             
               % computx x,y,z
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);

                
                IT=eye(3);
                xyz=[x;y;z];
                beta=-1i*delta^2*omega*sigma(mat(i))*muz/2;
                for iang=1:3
                    
                    uvec=IT(:,iang);
                    outcross=cross(uvec,xyz);
                    outcross2=cross(xyz,outcross);
                    
                    eiang=e(:,iang);
                    %coutcross=cross(xyz,eiang);
                    coutcross(1)=xyz(2)*eiang(3)-xyz(3)*eiang(2);
                    coutcross(2)=xyz(3)*eiang(1)-xyz(1)*eiang(3);
                    coutcross(3)=xyz(1)*eiang(2)-xyz(2)*eiang(1);
                    
                    for jj=1:3
                        cchkrank2(jj,iang)=cchkrank2(jj,iang)+det*intw(pp)*(coutcross(jj)+outcross2(jj))*beta/2;
                    end
                    
                    outcross=cross(uvec,xyz);
                    
                    ceiang=curle(:,iang);
                    
                    
                    for jj=1:3
                        nrank2(jj,iang)=nrank2(jj,iang)+det*intw(pp)*(uvec(jj)+ceiang(jj)/2)*(1-1/mu(mat(i)));
                    end
                    
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
                

                ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                        (ephy(pp,1:esizet)'*aeta(1:3))+...
                        (ephz(pp,1:esizet)'*azeta(1:3));
                
                % compute the solution for this problem, at this integration point
                % in this element
                e=zeros(3,3);
                
                e=ph'*lsol;

                
%                 % use stored functions
                ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                        (ecphy(pp,1:esizet)'*aseta(1:3))+...
                        (ecphz(pp,1:esizet)'*aszeta(1:3));
                
                % compute the solution for this problem, at this integration point
                % in this element
                curle=zeros(3,3);
                curle=ph'*lsol;
                
                               
               % computx x,y,z
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);

                IT=eye(3);
                xyz=[x;y;z];
                beta=-1i*delta^2*omega*sigma(mat(i))*muz/2;
                for iang=1:3
                    
                    uvec=IT(:,iang);
                    outcross=cross(uvec,xyz);
                    outcross2=cross(xyz,outcross);
                    
                    eiang=e(:,iang);
                    %coutcross=cross(xyz,eiang);
                    coutcross(1)=xyz(2)*eiang(3)-xyz(3)*eiang(2);
                    coutcross(2)=xyz(3)*eiang(1)-xyz(1)*eiang(3);
                    coutcross(3)=xyz(1)*eiang(2)-xyz(2)*eiang(1);
                    
                    
                    for jj=1:3
                        cchkrank2(jj,iang)=cchkrank2(jj,iang)+det*intw(pp)*(coutcross(jj)+outcross2(jj))*beta/2;
                    end
                    
                    outcross=cross(uvec,xyz);
                    
                    ceiang=curle(:,iang);
                    
                    
                    for jj=1:3
                        nrank2(jj,iang)=nrank2(jj,iang)+det*intw(pp)*(uvec(jj)+ceiang(jj)/2)*(1-1/mu(mat(i)));
                    end
                    
                    
                    
                end
            end
            
        end
        
        % end of loop over elements
    end
end
ptensor=-cchkrank2+nrank2;
format long e
ptensor=ptensor*delta^3

%exact=probdata.tensor


% ptensor=delta^3*(term1*eye(3)+term2)
% [av,ad]=eig(ptensor/delta^3)
% 
% % principal axis
% a11=1/sqrt(ad(1,1));
% a22=1/sqrt(ad(2,2));
% a33=1/sqrt(ad(3,3));
% 
% [x,y,z] = ellipsoid(0,0,0,a11,a22,a33);
% 
% [m n]=size(x);
% 
% % apply rotation matrix
% 
% for i=1:m
%     for j =1:n
%         xp=x(i,j);
%         yp=y(i,j);
%         zp=z(i,j);
%         
%         out=av'*[xp;yp;zp];
%         
%         xc(i,j)=out(1);
%         yc(i,j)=out(2);
%         zc(i,j)=out(3);
%     end
% end

% un comment for plot best fitting ellipsods
% mesh(xc,yc,zc);
% hidden off
% hold on
%
% [nt,dum]=size(myT)
%
% tetramesh(myT,coord,4*ones(nt,1))
% hold on


% ptensor_re=[   2.3220872576163054E-006  -1.0788772742365182E-007   7.9256122522951866E-008
%   -1.0788777587943266E-007   2.0459184841043421E-006   1.2858201502255533E-007
%    7.9256056334820725E-008   1.2858199854138293E-007   2.0955776835892287E-006];
% ptensor_im=[
%   -4.2234568236016352E-007  -1.1188450218761294E-008  -2.1652375590253432E-008
%   -1.1188127939360542E-008  -2.8648232341595222E-007  -1.0793700294431656E-008
%   -2.1652302007388236E-008  -1.0793689378120157E-008  -2.8031486930249233E-007];
% ptensor=-ptensor_im;
% [av,ad]=eig(ptensor/delta^3)
%
% % principal axis
% a11=1/sqrt(ad(1,1));
% a22=1/sqrt(ad(2,2));
% a33=1/sqrt(ad(3,3));
%
% [x,y,z] = ellipsoid(0,0,0,a11,a22,a33);
%
% [m n]=size(x);
%
% % apply rotation matrix
%
% for i=1:m
%     for j =1:n
%         xp=x(i,j);
%         yp=y(i,j);
%         zp=z(i,j);
%
%         out=av'*[xp;yp;zp];
%
%         xc(i,j)=out(1);
%         yc(i,j)=out(2);
%         zc(i,j)=out(3);
%     end
% end
%
% mesh(xc,yc,zc);
% hidden off
% hold on
%
%
% ptensor=ptensor_re;
% [av,ad]=eig(ptensor/delta^3)
%
% % principal axis
% a11=1/sqrt(ad(1,1));
% a22=1/sqrt(ad(2,2));
% a33=1/sqrt(ad(3,3));
%
% [x,y,z] = ellipsoid(0,0,0,a11,a22,a33);
%
% [m n]=size(x);
%
% % apply rotation matrix
%
% for i=1:m
%     for j =1:n
%         xp=x(i,j);
%         yp=y(i,j);
%         zp=z(i,j);
%
%         out=av'*[xp;yp;zp];
%
%         xc(i,j)=out(1);
%         yc(i,j)=out(2);
%         zc(i,j)=out(3);
%     end
% end
%
%
% mesh(xc,yc,zc);
% hidden off
% hold on
%
% [nt,dum]=size(myT)
%
% tetramesh(myT,coord,4*ones(nt,1))
% hold on
%
