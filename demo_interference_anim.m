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
    break;
end
    load(['meshes/h_' num2str(0.9^n) '_mesh.mat']);
    k = input('give k');
    tic;
    x = solver_interference(mesh,k);
    toc
    
    tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
    plot = trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x));
    xlabel('X'); ylabel('Y'); zlabel('real(u)');
    %set(plot,'FaceColor','none','EdgeColor','black');
    pause;
t=0;
while 1
    trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x.*exp(1i*t*k)));
    %set(plot,'ZData',real(x(tri))');%.*exp(1i*t*k)
    axis([-1 1 -1 1 -2 2]);
    drawnow;
    pause(0.2);
    t=t+0.05;
end