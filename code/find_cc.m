function [F1,V1,F2,V2] = find_cc(V,F)
%find_cc Finds the two connected components of the mesh defined by V and F
%   V are the vertices of the mesh, 3xN
%   F are the faces of the mesh, 3xM
% output is the faces and vertex of each component.
% This function assumes that there is only two connected components
    F = F';
    % Find the connected components
    % We must traverse all the vertices
    fSets = zeros(size(F,1),1,'uint32');
 
    % While some sets are not traversed
    cc = 0;
    while any(fSets==0)
        cc = cc + 1;
        fprintf('Connecting set #%d vertices...',cc);
        % Get next face available that is 0
        nextAvailFace = find(fSets==0,1,'first');
        % Get vertices of that face
        openVertices = F(nextAvailFace,:);
        while ~isempty(openVertices)
            % Find all indexs of the available faces
            availFaceInds = find(fSets==0);
            % Find available faces including any of the open vertices
            [availFaceSub, ~] = find(ismember(F(availFaceInds,:), openVertices));
            % Mark those faces as visited in fSets with the current set
            % number
            fSets(availFaceInds(availFaceSub)) = cc;
            % Open vertices are now the vertices of those marked faces.
            openVertices = F(availFaceInds(availFaceSub),:);
        end
        fprintf(' done! Set #%d has %d faces.\n',cc,nnz(fSets==cc));
    end
    
    % Get Vertices of cc 0
    setF = F(fSets==1,:);
    [unqVertIds, ~, newVertIndices] = unique(setF);
    F1 = reshape(newVertIndices,size(setF));
    V1 = V(:,unqVertIds);
    
    % Get Vertices of cc 1
    setF = F(fSets==2,:);
    [unqVertIds, ~, newVertIndices] = unique(setF);
    F2 = reshape(newVertIndices,size(setF));
    V2 = V(:,unqVertIds);    
end
