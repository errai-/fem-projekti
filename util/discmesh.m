% A function returning a mesh with maximum edge length h for 
% an origo-centered disc with radius r.

function mesh = discmesh(r, h)
    % Initialize the disc.
    disc = [1; 0; 0; r];
    % Create the geometry.
    g = decsg(disc);
    % Create the mesh and return.
    [p, e, t] = initmesh(g, 'Hmax', h,'MesherVersion', 'R2013a');
    mesh = inittri(p, t);
end

