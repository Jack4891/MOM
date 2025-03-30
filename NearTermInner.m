function I=NearTermInner(v1,v2,v3,rp,vFree,PlusOrMinus,TypeFlag)

%First find the free vertex
if v1==vFree; Vert1=v1; Vert2=v2; Vert3=v3; end
if v2==vFree; Vert1=v2; Vert2=v3; Vert3=v1; end
if v3==vFree; Vert1=v3; Vert2=v1; Vert3=v2; end

VertMat=[Vert1,Vert2,Vert3,Vert1];
r1=Vert1;
r2=Vert2;
r3=Vert3;
r12=r2-r1;
r13=r3-r1;
n=cross(r12,r13)./norm(cross(r12,r13));

%loop through all the observation points
for jj=1:length(rp(1,:))
    r=rp(:,jj);
    rho=r-n.*dot(n,r);
    d=dot(n,r-r1);
    rhoFree=vFree-n.*dot(n,vFree);
    
    for ii=1:3
        rhoP=VertMat(1:3,ii+1)-n.*(dot(n,VertMat(1:3,ii)));
        rhoM=VertMat(1:3,ii)-n.*(dot(n,VertMat(1:3,ii+1)));
        l=(rhoP-rhoM)./(norm(rhoP-rhoM));
        u=cross(l,n);
        lP=dot((rhoP-rho),l);
        lM=dot((rhoM-rho),l);
        P0=norm(dot((rhoP-rho),u));
        PP=norm(rhoP-rho);
        PM=norm(rhoM-rho);
        Psup0=((rhoP-rho)-lP.*l)./P0;
        Rsup0=sqrt(P0^2+d^2);
        RP=sqrt(PP^2+d^2);
        RM=sqrt(PM^2+d^2);
        if TypeFlag==1
            I1(1:3,ii)=((1/2)*(Rsup0^2*log((RP+lP)/(RM+lM))+lP*RP-lM*RM)).*(u);
            I2(1:3,ii)=(dot(Psup0,u)*(P0*log((RP+lP)/(RM+lM))-abs(d)*...
            (atan2(P0*lP,Rsup0^2+abs(d)*RP)-atan2(P0*lM,Rsup0^2+abs(d)*RM)))).*(rho-rhoFree);
        
        elseif TypeFlag==2
            I2(ii)=(dot(Psup0,u)*(P0*log((RP+lP)/(RM+lM))-abs(d)*...
            (atan2(P0*lP,Rsup0^2+abs(d)*RP)-atan2(P0*lM,Rsup0^2+abs(d)*RM))));
        end
    end
    if TypeFlag==1
        if PlusOrMinus==2
            I(1:3,jj)=sum(I1,2)+sum(I2,2);
        else
            I(1:3,jj)=-sum(I1,2)-sum(I2,2);
        end
    elseif TypeFlag==2
        I(1,jj)=sum(I2);
    end
end