clc
clear
close all

freq = 10 * 1e9;
Lambda = 3e8/freq;
w0 = 2*pi*freq;

mu0 = 4*pi*1e-7;
e0 = 8.854e-12;

k0 = w0 * sqrt(mu0*e0);
L = Lambda;
W = Lambda;
% L = 2;
% W = 3;


%% generate the mesh
n = 10; % number of elements along the x axis
x = linspace(-L/2,L/2,n);
y = linspace(-W/2,W/2,n);

[X,Y] = meshgrid(x,y);

X = X(:);
Y = Y(:);

DT = delaunay(X,Y);

figure
triplot(DT,X,Y)
xlabel("x-axis")
ylabel("y-axis")
title("triangular mesh")

Edges = [];
edgesTotal = 0;
for n = 1:1:size(DT,1) 
    % disp(DT(n,:)) 
    % DT(n,:) refers to the triangle vertexes in row n
    tri1 = DT(n,1:3)';
    for k = (n + 1):1:size(DT,1)
        tri2 = DT(k,1:3)';
        a=1-all([tri1-tri2(1) tri1-tri2(2) tri1-tri2(3)]); 
        if(sum(a) == 2)
            edgesTotal = edgesTotal + 1;
            Edges = [Edges tri2((find(a)))];
            posTri(edgesTotal) = n;  % using edgesTotal as an index
            negTri(edgesTotal) = k;
        end
    end
end

numberTriangles = size(DT,1);
Triangles = DT';  

p = [X Y zeros(length(X),1)]';

for m = 1:numberTriangles
    N = Triangles(1:3,m);
    Vec1 = p(:,N(1)) - p(:,N(2));
    Vec2 = p(:,N(3)) - p(:,N(2));

    Area(m) = norm(cross(Vec1,Vec2))/2;

    Center(:,m) = 1/3 * sum(p(:,N),2);
end

for m = 1:edgesTotal
    edgeLength(m) = norm(p(:,Edges(1,m)) - p(:,Edges(2,m)));
end

nearDist = mean(edgeLength);

numEdges = size(Edges(1,:));

Triangles(4,:) = 1;

% Plus triangle position
for m = 1:edgesTotal
    triIndex = posTri(m);
    n1 = Triangles(1,triIndex);
    n2 = Triangles(2,triIndex);
    n3 = Triangles(3,triIndex);
    if((n1 ~= Edges(1,m) &&(n1 ~= Edges(2,m)))) NODE = n1;end
    if((n2 ~= Edges(1,m) &&(n2 ~= Edges(2,m)))) NODE = n2;end
    if((n3 ~= Edges(1,m) &&(n3 ~= Edges(2,m)))) NODE = n3;end
    freeVertexPos(:,m) = p(:,NODE);
    rhoPlus(:,m) = -Center(:,triIndex) + freeVertexPos(:,m);
end

% negative triangle position
for m = 1:edgesTotal
    triIndex = negTri(m);
    n1 = Triangles(1,triIndex);
    n2 = Triangles(2,triIndex);
    n3 = Triangles(3,triIndex);
    if((n1 ~= Edges(1,m) &&(n1 ~= Edges(2,m)))) NODE = n1;end
    if((n2 ~= Edges(1,m) &&(n2 ~= Edges(2,m)))) NODE = n2;end
    if((n3 ~= Edges(1,m) &&(n3 ~= Edges(2,m)))) NODE = n3;end
    freeVertexNeg(:,m) = p(:,NODE);
    rhoMinus(:,m) = Center(:,triIndex) - freeVertexNeg(:,m);
end

%% build the impedance matrix
numQRuleO = 1;
numQRuleS = 1;

QRuleS = GetQuadRule(numQRuleO);
QRuleO = GetQuadRule(numQRuleO);

for m = 1:size(Edges,2)
    % m+ 
    V1mp = p(:,Triangles(1,posTri(m)));
    V2mp = p(:,Triangles(2,posTri(m)));
    V3mp = p(:,Triangles(3,posTri(m)));
    % m+ 
    V1mm = p(:,Triangles(1,negTri(m)));
    V2mm = p(:,Triangles(2,negTri(m)));
    V3mm = p(:,Triangles(3,negTri(m)));

    % rp and rhom vectors
    rpm = repmat(QRuleO(:,1),[1,3]).'.*repmat(V1mm,[1,numQRuleO])+...
        repmat(QRuleO(:,2),[1,3]).'.*repmat(V2mm,[1,numQRuleO])+...
        repmat(QRuleO(:,3),[1,3]).'.*repmat(V3mm,[1,numQRuleO]);
    RHOm_Minus = rpm - repmat(freeVertexNeg(:,m),[1,numQRuleO]); 

    rpp = repmat(QRuleO(:,1),[1,3]).'.*repmat(V1mp,[1,numQRuleO])+...
    repmat(QRuleO(:,2),[1,3]).'.*repmat(V2mp,[1,numQRuleO])+...
    repmat(QRuleO(:,3),[1,3]).'.*repmat(V3mp,[1,numQRuleO]);
    RHOm_Plus = repmat(freeVertexPos(:,m),[1,numQRuleO]) - rpp;

    for n = 1:size(Edges,2)
        % triangle n+
        V1np = p(:,Triangles(1,posTri(n)));
        V2np = p(:,Triangles(2,posTri(n)));
        V3np = p(:,Triangles(3,posTri(n)));
        % triangle n-
        V1nm = p(:,Triangles(1,negTri(n)));
        V2nm = p(:,Triangles(2,negTri(n)));
        V3nm = p(:,Triangles(3,negTri(n)));

        % rq vectors and rhon vectors
        rqm = repmat(QRuleS(:,1),[1,3]).'.*repmat(V1nm,[1,numQRuleS])+...
            repmat(QRuleS(:,2),[1,3]).'.*repmat(V2nm,[1,numQRuleS])+...
            repmat(QRuleS(:,3),[1,3]).'.*repmat(V3nm,[1,numQRuleS]);
        RHOn_Minus = rqm - repmat(freeVertexNeg(:,n),[1,numQRuleS]); 
        
        RHO_Saved_Minus(:,n) = RHOn_Minus;
    
        rqp = repmat(QRuleS(:,1),[1,3]).'.*repmat(V1np,[1,numQRuleS])+...
            repmat(QRuleS(:,2),[1,3]).'.*repmat(V2np,[1,numQRuleS])+...
            repmat(QRuleS(:,3),[1,3]).'.*repmat(V3np,[1,numQRuleS]);
        RHOn_Plus = repmat(freeVertexPos(:,n),[1,numQRuleS]) - rqp; 

        RHO_Saved_Plus(:,n) = RHOn_Plus;


        %% Tm- and Tn-
        centroid_Tmm = 1/3 * sum(p(:,Triangles(1:3,negTri(m))),2);
        centroid_Tnm = 1/3 * sum(p(:,Triangles(1:3,negTri(n))),2);
        deltaCentroid = norm(centroid_Tnm-centroid_Tmm);
        % not a self term
        Rpq = sqrt(pdist2(rpm(1,:).',rqm(1,:).').^2 + pdist2(rpm(2,:).',rqm(2,:).').^2 + ...
                pdist2(rpm(3,:).',rqm(3,:).').^2 );
        if(deltaCentroid>nearDist)

            TmmTnm=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOn_Minus.'*RHOm_Minus-...
            (1i./(w0.*e0))).*(exp(-1i.*k0.*Rpq)./Rpq))));

        else
            % term 1
            Term1=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOm_Minus.'*RHOn_Minus-(1i./(w0.*e0))).*...
            GFuncNonSingPart(Rpq,k0))));
            % term 2
            Term2=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(negTri(n)))).*...
            QRuleO(:,4).'.*((1i.*w0.*mu0./4).*...
            dot(RHOm_Minus,NearTermInner(V1nm,V2nm,V3nm,rpm,freeVertexNeg(:,n),2,1))));
            % term 3
            Term3=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(negTri(n)))).*...
            QRuleO(:,4).'.*((-1i./(w0.*e0)).*NearTermInner(V1nm,V2nm,V3nm,rpm,freeVertexNeg(:,n),2,2)));
            % add 
            TmmTnm=Term1+Term2+Term3; 
        end

        %% Tm- and Tn+
        centroid_Tmm = 1/3 * sum(p(:,Triangles(1:3,negTri(m))),2);
        centroid_Tnp = 1/3 * sum(p(:,Triangles(1:3,posTri(n))),2);
        deltaCentroid = norm(centroid_Tnp-centroid_Tmm);
        % not a self term
        Rpq = sqrt(pdist2(rpm(1,:).',rqp(1,:).').^2 + pdist2(rpm(2,:).',rqp(2,:).').^2 + ...
                pdist2(rpm(3,:).',rqp(3,:).').^2 );
        if(deltaCentroid>nearDist)

            TmmTnp=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOn_Plus.'*RHOm_Minus + ...
            (1i./(w0.*e0))).*(exp(-1i.*k0.*Rpq)./Rpq))));

        else
            % term 1
            Term1=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOm_Minus.'*RHOn_Plus + (1i./(w0.*e0))).*...
            GFuncNonSingPart(Rpq,k0))));
            % term 2
            Term2=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(posTri(n)))).*...
            QRuleO(:,4).'.*((1i.*w0.*mu0./4).*...
            dot(RHOm_Minus,NearTermInner(V1np,V2np,V3np,rpm,freeVertexPos(:,n),2,1))));
            % term 3
            Term3=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(posTri(n)))).*...
            QRuleO(:,4).'.*((1i./(w0.*e0)).*NearTermInner(V1np,V2np,V3np,rpm,freeVertexPos(:,n),2,2)));
            % add 
            TmmTnp=Term1+Term2+Term3; 
        end

        %% Tm+ and Tn-
        centroid_Tmp = 1/3 * sum(p(:,Triangles(1:3,posTri(m))),2);
        centroid_Tnm = 1/3 * sum(p(:,Triangles(1:3,negTri(n))),2);
        deltaCentroid = norm(centroid_Tnm-centroid_Tmp);
        % not a self term
        Rpq = sqrt(pdist2(rpp(1,:).',rqm(1,:).').^2 + pdist2(rpp(2,:).',rqm(2,:).').^2 + ...
                pdist2(rpp(3,:).',rqm(3,:).').^2 );
        if(deltaCentroid>nearDist)

            TmpTnm=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOn_Minus.'*RHOm_Plus + ...
            (1i./(w0.*e0))).*(exp(-1i.*k0.*Rpq)./Rpq))));

        else
            % term 1
            Term1=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOm_Plus.'*RHOn_Minus + (1i./(w0.*e0))).*...
            GFuncNonSingPart(Rpq,k0))));
            % term 2
            Term2=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(negTri(n)))).*...
            QRuleO(:,4).'.*((1i.*w0.*mu0./4).*...
            dot(RHOm_Plus,NearTermInner(V1nm,V2nm,V3nm,rpp,freeVertexNeg(:,n),2,1))));
            % term 3
            Term3=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(negTri(n)))).*...
            QRuleO(:,4).'.*((1i./(w0.*e0)).*NearTermInner(V1nm,V2nm,V3nm,rpp,freeVertexNeg(:,n),2,2)));
            % add 
            TmpTnm=Term1+Term2+Term3; 
        end

        %% Tm+ and Tn+
        centroid_Tmp = 1/3 * sum(p(:,Triangles(1:3,posTri(m))),2);
        centroid_Tnp = 1/3 * sum(p(:,Triangles(1:3,posTri(n))),2);
        deltaCentroid = norm(centroid_Tnp-centroid_Tmp);
        % not a self term
        Rpq = sqrt(pdist2(rpp(1,:).',rqp(1,:).').^2 + pdist2(rpp(2,:).',rqp(2,:).').^2 + ...
                pdist2(rpp(3,:).',rqp(3,:).').^2 );
        if(deltaCentroid>nearDist)

            TmpTnp=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOn_Plus.'*RHOm_Plus-...
            (1i./(w0.*e0))).*(exp(-1i.*k0.*Rpq)./Rpq))));

        else
            % term 1
            Term1=sum(sum(((edgeLength(m).*edgeLength(n)./(4.*pi)).*...
            (QRuleO(:,4)*QRuleS(:,4).').*((1i.*w0.*mu0./4).*RHOm_Plus.'*RHOn_Plus-(1i./(w0.*e0))).*...
            GFuncNonSingPart(Rpq,k0))));
            % term 2
            Term2=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(posTri(n)))).*...
            QRuleO(:,4).'.*((1i.*w0.*mu0./4).*...
            dot(RHOm_Plus,NearTermInner(V1np,V2np,V3np,rpp,freeVertexPos(:,n),2,1))));
            % term 3
            Term3=sum((edgeLength(m).*edgeLength(n)./(4.*pi.*Area(posTri(n)))).*...
            QRuleO(:,4).'.*((-1i./(w0.*e0)).*NearTermInner(V1np,V2np,V3np,rpp,freeVertexPos(:,n),2,2)));
            % add 
            TmpTnp=Term1+Term2+Term3; 
        end


        %% ZMN
        Zmn(m,n) = TmmTnm + TmmTnp + TmpTnp + TmpTnm;
    end     


    %% generation of Vm
    IncTheta=90*pi/180;
    IncPhi=0*pi/180;
    IncK=[sin(IncTheta)*cos(IncPhi) sin(IncTheta)*sin(IncPhi) ... 
    cos(IncTheta)]; 
    IncPol=[1;0;0];
    E0i=1; 
    Einc=@(x,y,z)E0i.*IncPol.*exp(-1i.*dot(k0.*IncK,[x,y,z]));


    V(m,1)=sum((edgeLength(m)./2).*QRuleO(:,4).'.*dot(RHOm_Minus,Einc(rpm(1,:),rpm(2,:),rpm(3,:))))+...
        sum((edgeLength(m)./2).*QRuleO(:,4).'.*dot(RHOm_Plus,Einc(rpp(1,:),rpp(2,:),rpp(3,:))));
end

%% currents
I = inv(Zmn) * V

%% plot surface current densities
Index=find(Triangles(4,:)<=1);
numTriangles=length(Index);

for k = 1:numTriangles
    J = [0 0 0];
    for m = 1:edgesTotal
        IE = I(m) * edgeLength(m);
        if(posTri(m) == k)
            J = J + IE * RHO_Saved_Plus(:,m)' / (2*Area(posTri(m)));
        end
        if(negTri(m) == k)
            J = J + IE * RHO_Saved_Minus(:,m)' / (2*Area(posTri(m)));
        end
        Jx = J(1);
        Jy = J(2);
        Jz = J(3);
        normalCurrent(k) = abs(norm(J));


    end
end

Jmax=max(normalCurrent);
MaxCurrent=strcat(num2str(Jmax),'[A/m]')
CurrentNorm1=normalCurrent./max(normalCurrent);
for m=1:numTriangles
    N=Triangles(1:3,m);
    Xp(1:3,m)=[p(1,N)]';
    Yp(1:3,m)=[p(2,N)]';
    Zp(1:3,m)=[p(3,N)]';
end
Cam=repmat(CurrentNorm1,3,1);
Cph=repmat(Jy./max(normalCurrent),3,1);
figure
h = fill3(Xp,Yp,Zp,Cam);
colormap jet;
axis('equal')

xlim([-L L]) 
ylim([-W W])
