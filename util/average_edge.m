function h = longest_edge( mesh )
%LONGEST_EDGE Finds the length of the longest edge

        px = mesh.p(1,:);    % 1xNpoint - vector
        py = mesh.p(2,:);    % 1xNpoint - vector
        ex = px(mesh.edges); % 2xNedges - matrix.
        ey = py(mesh.edges); % 2xNedges - matrix.
        
        % pythagoras
        len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;
        h = mean(sqrt(len));
        
end

