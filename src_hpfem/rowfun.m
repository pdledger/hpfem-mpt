function out=rowfun(unkz,unksid,unkint,order,glob,globfa,...
                        i,unkfatp1,unkfatp2,unkfatp3,esize)
                    
                    
nbas=0;
out = zeros(esize,1);
% edge basis functions
for m=1:6
    nbas=nbas+1;
    out(nbas)=unkz(glob(i,m));
end

if order>=1
    for pp=1:order
        for m=1:6
            nbas=nbas+1;
            out(nbas)=unksid(glob(i,m),pp);
        end
    end
end

% face basis functions
if order>=2
    
    for m=1:4

        for k=1:(order*order-order)/2
              nbas=nbas+1;
              out(nbas)=unkfatp1(globfa(i,m),k);
        end
    end
    
    for m=1:4

        for k=1:(order*order-order)/2
              nbas=nbas+1;
              out(nbas)=unkfatp2(globfa(i,m),k);
        end
    end
    for m=1:4
        for jj=0:order-2
            nbas=nbas+1;
            out(nbas)=unkfatp3(globfa(i,m),jj+1);
        end
    end

    end
if order>= 3
% interior functions
    k=0;
    for ii=0:order-3
        for jj=0:order-3
            for kk=0:order-3
                if ii+jj+kk <= order-3
                    nbas=nbas+1;
                    k=k+1;
                    out(nbas)=unkint(i,k);
                end
            end
        end
    end
    
    %     non-gradients
    nintbas=k;
    for ii=nintbas+1:(order-2)*(order-1)*(order+1)/2
        k=k+1;
        nbas=nbas+1;
        out(nbas)=unkint(i,k);
        
    end
%     for ii=0:order-3
%         for jj=0:order-3
%             for kk=0:order-3
%                 if ii+jj+kk<=order-3
%                     k=k+1;
%                     %     typ 2a
%                     nbas=nbas+1;
%                     out(nbas)=unkint(i,k);
%                     k=k+1;
%                     %     typ 2b
%                     nbas=nbas+1;
%                     out(nbas)=unkint(i,k);
%                 end
%             end
%         end
%     end
%     
%     for jj=0:order-3
%         for kk=0:order-3
%             if jj+kk<=order-3
%                 k=k+1
%                 nbas=nbas+1
%                 out(nbas)=unkint(i,k);
%             end
%         end
%     end
    
    
    
    %     for j=1:(order-2)*(order-1)*(order+1)/2
    %         nbas=nbas+1;
    %         out(nbas)=unkint(i,j);
    %     end
end