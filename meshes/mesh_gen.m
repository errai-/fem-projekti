time=0;
while(time<50)
    tic;
    mesh=discmesh(1,h);
    time=toc
    save(['h_',num2str(h),'_mesh.mat'], 'mesh');
    h=0.9*h
end
