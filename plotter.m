clear all;
load('data.mat');

xdata=[];
for datum=data
    xdata=[xdata repmat(datum.h_avg, size(datum.k))];
end
ydata=[data.k];
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('h (average)'); ylabel('k'); zlabel('total error');

pause;

xdata=[];
for datum=data
    xdata=[xdata repmat(datum.h_max, size(datum.k))];
end
ydata=[data.k];
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('h (max)'); ylabel('k'); zlabel('total error');

pause;

xdata=[];
for datum=data
    xdata=[xdata repmat(datum.numnodes, size(datum.k))];
end
ydata=[data.k];
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('number of nodes'); ylabel('k'); zlabel('total error');

pause;

xdata=[];
for datum=data
    xdata=[xdata repmat(datum.numnodes, size(datum.k))];
end
ydata=[data.k].^3;
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('number of nodes'); ylabel('k^3'); zlabel('total error');

pause;
xdata=[];
for datum=data
    xdata=[xdata datum.h_avg*datum.k];
end
ydata=[data.k];
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('h (average) *k'); ylabel('k'); zlabel('total error');

pause;

xdata=[];
for datum=data
    xdata=[xdata datum.h_avg*(datum.k.^1.5)];
end
ydata=[data.k];
zdata=[data.error];
tri = delaunay(xdata, ydata);
trisurf(tri, xdata, ydata, zdata);
set(gca,'xscale','log','yscale','log','zscale','log');
xlabel('h (average) *k^1.5'); ylabel('k'); zlabel('total error');
