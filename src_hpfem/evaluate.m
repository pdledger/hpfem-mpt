function [ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,ecphy1,ecphz1,ecphx2,...
          ecphy2,ecphz2,gphx1,gphy1,gphz1,gphx2,gphy2,gphz2,phh11,phh12]= evaluate(nip,intxi,inteta,intzeta,order,esizet,gorder)

ephx1 = zeros(nip,esizet);
ephy1 = zeros(nip,esizet);
ephz1 = zeros(nip,esizet);

ephx2 = zeros(nip,esizet);
ephy2 = zeros(nip,esizet);
ephz2 = zeros(nip,esizet);

ecphx1 = zeros(nip,esizet);
ecphy1 = zeros(nip,esizet);
ecphz1 = zeros(nip,esizet);

ecphx2 = zeros(nip,esizet);
ecphy2 = zeros(nip,esizet);
ecphz2 = zeros(nip,esizet);

gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;

gphx1=zeros(nip,gesizet);
gphy1=zeros(nip,gesizet);
gphz1=zeros(nip,gesizet);

gphx2=zeros(nip,gesizet);
gphy2=zeros(nip,gesizet);
gphz2=zeros(nip,gesizet);

phh11=zeros(nip,gesizet);
phh12=zeros(nip,gesizet);


for i=1:nip
    % set up mapping functions
    axi(1)=1;
    axi(2)=0;
    axi(3)=0;
    
    aeta(1)=0;
    aeta(2)=1;
    aeta(3)=0;
    
    azeta(1)=0;
    azeta(2)=0;
    azeta(3)=1;
    
    asxi(1)=1;
    asxi(2)=0;
    asxi(3)=0;
    
    aseta(1)=0;
    aseta(2)=1;
    aseta(3)=0;
    
    aszeta(1)=0;
    aszeta(2)=0;
    aszeta(3)=1;
    % write(6,*)intxi(i),inteta(i),intzeta(i),i,nip
    ph= basis(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,order,1,esizet);
    
    % write(6,*)'evaluated basis',i
    ephx1(i,1:esizet)=ph(1:esizet,1)';
    ephy1(i,1:esizet)=ph(1:esizet,2)';
    ephz1(i,1:esizet)=ph(1:esizet,3)';
    
    ph= basis(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,order,2,esizet);
    
    % write(6,*)'evaluated basis typ 2',i
    ephx2(i,1:esizet)=ph(1:esizet,1)';
    ephy2(i,1:esizet)=ph(1:esizet,2)';
    ephz2(i,1:esizet)=ph(1:esizet,3)';
    
    ph= curlbasis(intxi(i),inteta(i),intzeta(i),asxi,aseta,aszeta,order,1,esizet);
    ecphx1(i,1:esizet)=ph(1:esizet,1)';
    ecphy1(i,1:esizet)=ph(1:esizet,2)';
    ecphz1(i,1:esizet)=ph(1:esizet,3)';
    
    ph= curlbasis(intxi(i),inteta(i),intzeta(i),asxi,aseta,aszeta,order,2,esizet);
    ecphx2(i,1:esizet)=ph(1:esizet,1)';
    ecphy2(i,1:esizet)=ph(1:esizet,2)';
    ecphz2(i,1:esizet)=ph(1:esizet,3)';
    
    
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,1,gesizet);

    gphx1(i,1:gesizet)=gph(1:gesizet,1)';
    gphy1(i,1:gesizet)=gph(1:gesizet,2)';
    gphz1(i,1:gesizet)=gph(1:gesizet,3)';
    
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,2,gesizet);
    gphx2(i,1:gesizet)=gph(1:gesizet,1)';
    gphy2(i,1:gesizet)=gph(1:gesizet,2)';
    gphz2(i,1:gesizet)=gph(1:gesizet,3)';

    
    phh1=basish1(gesizet,intxi(i),inteta(i),intzeta(i),gorder+1,1);
    phh11(i,1:gesizet)=phh1(1:gesizet,1)';
    
    phh1=basish1(gesizet,intxi(i),inteta(i),intzeta(i),gorder+1,2);
    phh12(i,1:gesizet)=phh1(1:gesizet,1)';
    
    
    %     l glagbasis(phlag,lagesizet,intxi(i),inteta(i),intzeta(i),
%     !     &axi,aeta,azeta,1,order)
%     !      do j=1,lagesizet
%     !      lphx1(i,j)=phlag(j,1)
%     !      lphy1(i,j)=phlag(j,2)
%     !      lphz1(i,j)=phlag(j,3)
%     !      enddo
%     !
%     !      call glagbasis(phlag,lagesizet,intxi(i),inteta(i),intzeta(i),
%     !     &axi,aeta,azeta,2,order)
%     !      do j=1,lagesizet
%     !      lphx2(i,j)=phlag(j,1)
%     !      lphy2(i,j)=phlag(j,2)
%     !      lphz2(i,j)=phlag(j,3)
%     !      enddo 
    
end