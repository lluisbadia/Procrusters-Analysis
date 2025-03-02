%% Left pair registration
% Load the .vtk files from two subjects
subject1 = 'path\to\file1';
subject2 = 'path\to\file2';
[left_vertex1, left_face1] = read_vtk(subject1, 1);
[left_vertex2, left_face2] = read_vtk(subject2, 1);

VA = left_vertex1'; % Convert from 3 x N to N x 3
VB = left_vertex2'; % Convert from 3 x N to N x 3

% Perform Procrustes registration
% Align subject 2 left hippocampus to subject 1 left hippocampus
[d, Z] = my_procrustes(VA, VB);
fprintf('Procrustes distance: %f\n', d);

% Visualize the result
redColors = repmat([1, 0, 0], size(VA,1), 1); % Red for Subject 1
blueColors = repmat([0, 0, 1], size(VB,1), 1); % Blue for Subject 2
figure;
plot_mesh(VA, left_face1, struct('face_vertex_color', redColors));
hold on;
plot_mesh(VB, left_face1, struct('face_vertex_color', blueColors));
title('Left Hippocampus Meshes Before Alignment');

redColors = repmat([1, 0, 0], size(VA,1), 1); % Red for Subject 1
blueColors = repmat([0, 0, 1], size(Z,1), 1); % Blue for Subject 2
figure;
plot_mesh(VA, left_face1, struct('face_vertex_color', redColors));
hold on;
plot_mesh(Z, left_face1, struct('face_vertex_color', blueColors));
title('Left Hippocampus Meshes After Alignment');

%% Right pair registration
% Load the .vtk files from two subjects
subject1 = 'path\to\file1';
subject2 = 'path\to\file2';
[right_vertex1, right_face1] = read_vtk(subject1, 1);
[right_vertex2, right_face2] = read_vtk(subject2, 1);

VA = right_vertex1'; % Convert from 3 x N to N x 3
VB = right_vertex2'; % Convert from 3 x N to N x 3

% Perform Procrustes registration
% Align subject 2 right hippocampus to subject 1 right hippocampus
[d, Z] = my_procrustes(VA, VB);
fprintf('Procrustes distance: %f\n', d);

% Visualize the result
redColors = repmat([1, 0, 0], size(VA,1), 1); % Red for Subject 1
blueColors = repmat([0, 0, 1], size(VB,1), 1); % Blue for Subject 2
figure;
plot_mesh(VA, right_face1, struct('face_vertex_color', redColors));
hold on;
plot_mesh(VB, right_face1, struct('face_vertex_color', blueColors));
title('Right Hippocampus Meshes Before Alignment');

redColors = repmat([1, 0, 0], size(VA,1), 1); % Red for Subject 1
blueColors = repmat([0, 0, 1], size(Z,1), 1); % Blue for Subject 2
figure;
plot_mesh(VA, right_face1, struct('face_vertex_color', redColors));
hold on;
plot_mesh(Z, right_face1, struct('face_vertex_color', blueColors));
title('Right Hippocampus Meshes After Alignment');