% A function returning a mesh with maximum edge length h for 
% an origo-centered disc with radius r.

function mesh = discmesh(r, h)
    % Initialize the disc.
    disc = [1; 0; 0; r];
    % Create the geometry.
    g = descg(disc);
    % Create the mesh and return.
    [p, e, t] = initmesh(g, 'Hmax', h);
    mesh = inittri(p, t)
end
