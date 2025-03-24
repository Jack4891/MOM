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
    tri1 = DT(n,:);
    for k = (n + 1):1:size(DT,1)
        tri2 = DT(k,:);
        a=[tri1-tri2(1) tri1-tri2(2) tri1-tri2(3)] == zeros(1,9);
        if(sum(a) == 2)
            edgesTotal = edgesTotal + 1;
            Edges = [Edges tri1((find(a)))];
            posTri(edgesTotal) = Tri1;
            negTri(edgesTotal) = Tri2;

        end
    end
end

