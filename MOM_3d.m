clc
clear
close all

freq = 10 * 1e9;
Lambda = 3e8/freq;

% L = 20 * Lambda;
% W = 10 * Lambda;
L = 2;
W = 3;

n = 3; % number of elements along the x axis
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
    disp("your mom")
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



