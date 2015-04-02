clear all;
addpath('./util');

mesh_files = dir('meshes/*.mat'); 
numfiles = length(mesh_files);
while 1
    n = input('give mesh index');
    if(exist(['meshes/h_' num2str(0.9^n) '_mesh.mat'], 'file')==0)
        disp('no mesh found');
        continue;
    end
    load(['meshes/h_' num2str(0.9^n) '_mesh.mat']);
    k = input('give k');
    tic;
    [x,err] = solver(mesh,k);
    toc
    
    tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
    trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x));
    xlabel('X'); ylabel('Y'); zlabel('real(u)');
    pause;
    
    tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
    trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', abs(x));
    xlabel('X'); ylabel('Y'); zlabel('abs(u)');

end