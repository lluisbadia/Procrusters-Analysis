% Specify folder containing VTK files:
meshFolder = 'path\to\folder';
% Get all files from the meshFolder matching out_mesh_*.vtk:
files = dir(fullfile(meshFolder, 'out_mesh_*.vtk'));

% Loop through each file
for i = 1:length(files)
    % Full path to the current file
    fileName = fullfile(meshFolder, files(i).name);
    fprintf('Processing file: %s\n', fileName);

    % Read the mesh with read_vtk function
    % 'vertex' is a 'N vertices x 3' array specifying vertex position
    % 'face' is a 'M faces x 3' array specifying mesh connectivity
    [vertex, face] = read_vtk(fileName);

    % Check that vertices and faces have the correct size
    [vertex, face] = check_face_vertex(vertex, face);

    % Separate both meshes into connected components
    % V are the vertices of the mesh, 3xN
    % F are the faces of the mesh, 3xM
    [F1, V1, F2, V2] = find_cc(vertex, face);

    % Determine if it is left or right
    % Anatomically, left side is located more to the left, lower mean
    meanX1 = mean(V1(1,:));
    meanX2 = mean(V2(1,:));
    % Assign left or right side given the computed x-coordinate means
    if meanX1 < meanX2
        left_vertex = V1;
        left_face = F1;
        right_vertex = V2;
        right_face = F2;
    else
        left_vertex = V2;
        left_face = F2;
        right_vertex = V1;
        right_face = F1;
    end

    % Save as VTK using the write_vtk function
    [~, baseName] = fileparts(files(i).name);
    leftVTKName  = fullfile(meshFolder, [baseName '_left.vtk']);
    rightVTKName = fullfile(meshFolder, [baseName '_right.vtk']);
    write_vtk(left_vertex, left_face, leftVTKName, 1);
    write_vtk(right_vertex, right_face, rightVTKName, 1);
    fprintf('  -> Saved left and right hippocampus of %s\n', baseName);
end