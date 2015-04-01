clear all;
addpath('./util');

mesh_files = dir('meshes/*.mat'); 
numfiles = length(mesh_files);
load('data.mat');

for mesh_idx = 1:numfiles 
    load(['meshes/' mesh_files(mesh_idx).name]);
    datum = struct('h_max',longest_edge(mesh),'h_avg',average_edge(mesh),'numnodes',size(mesh.p,2),'error',[],'k',[]);
    datum.k = 1*(1.1).^(61:1:90);
    for k = datum.k
        tic;
        [x,err] = solver(mesh,k);
        time=toc;
        datum.error = [datum.error err];
        disp(time);
    end
    data=[data datum]
end
save('data.mat', 'data');
