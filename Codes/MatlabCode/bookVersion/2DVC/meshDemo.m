% Mesh generation, refinement and visualization

hold off; close all; clear all;

%% Create a L-shaped domain

% Specifying the coordinates.
nodeL = [0,-1;1,-1;1,0;1,1;0,1;-1,1;-1,0;0,0];  

% Specifying the connectivity.
elemL = [1,2,8;3,8,2;8,3,5;4,5,3;7,8,6;5,6,8]; 
figure; hold off
showmesh(nodeL,elemL); 
axis on;

% Plotting the indices of all triangles.
findelem(nodeL,elemL); 

% Plotting the indices of all vertices.
findnode(nodeL);       
title('Mesh of L-shaped domain',...
  'Units','normalized','Position',[0.5,1.05,0]);


%% Create a rectangular domain with mesh size specified

% Rectangle [0,1]*[0,2], mesh size 0.5
[nodeRectangular,elemRectangular] = squaremesh([0,1,0,2],0.5);
figure; hold off
showmesh(nodeRectangular,elemRectangular); 
axis on;
findelem(nodeRectangular,elemRectangular);  
findnode(nodeRectangular);      
title('Mesh of rectangular domain',...
  'Units','normalized','Position',[0.5,1.05,0]);

%% Create a circular domain with mesh size specified

% Circular mesh with center (0,1) and radius 2, mesh size 0.5
[nodeRectangular,elemRectangular] = circlemesh(0,1,2,0.5);
figure; hold off
showmesh(nodeRectangular,elemRectangular); 
axis on;
findelem(nodeRectangular,elemRectangular);  
findnode(nodeRectangular);      
title('Mesh of circular domain',...
  'Units','normalized','Position',[0.5,1.05,0]);

%% Create a customized domain

% Cat domain created with Im2Mesh package
% Jiexian Ma (2025). Im2mesh (2D image to triangular meshes) 
% (https://www.mathworks.com/matlabcentral/fileexchange/71772
% -im2mesh-2d-image-to-triangular-meshes),
% MATLAB Central File Exchange. 
% Retrieved January 30, 2025.

load("cat.mat")
figure; hold off
showmesh(nodeCat,elemCat); 
axis on;      
title('Mesh of cat domain',...
  'Units','normalized','Position',[0.5,1.05,0]);

%% Uniform refinement of meshes

% Refine mesh uniformly by quadrisecting
[nodeRefinedL,elemRefinedL] = uniformrefine(nodeL,elemL); 

figure; hold off
showmesh(nodeRefinedL,elemRefinedL); 
axis on;
findelem(nodeRefinedL,elemRefinedL); 
findnode(nodeRefinedL);       
title('Mesh of refined L-shaped domain',...
  'Units','normalized','Position',[0.5,1.05,0]);
hold off