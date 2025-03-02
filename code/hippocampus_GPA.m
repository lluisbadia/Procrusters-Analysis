%% Process Left Hippocampus
% Select one shape to be the initial approximate mean shape
folder = 'path\to\folder';
files = dir(fullfile(folder, '*_left.vtk')); % Get files with _left.vtk
numShapes = numel(files);
shapes = cell(1, numShapes);

% Load each .vtk file and extract the vertex data (transposed to Nx3)
disp('Loading left hippocampus files...');
for i = 1:numShapes
    fileName = fullfile(folder, files(i).name);
    [vertex, face] = read_vtk(fileName, 0); % Read .vtk file
    shapes{i} = vertex'; % Transpose
end

% Set the initial approximate mean shape
initial_shape = shapes{1};
reference = initial_shape;

% Set convergence parameters for GPA
maxIter = 10; % Maximum number of iterations
iter = 0;
threshold = 1e-10;

% GPA Iteration for Left Hippocampus:
while iter < maxIter
    iter = iter + 1;
    alignedShapes = cell(1, numShapes);

    % Align each shape to the current reference using Procrustes
    for i = 1:numShapes
        [~, Z] = my_procrustes(reference, shapes{i});
        alignedShapes{i} = Z;
    end

    % Calculate the new approximate mean shape from the aligned shapes
    meanShape = zeros(size(reference));
    for i = 1:numShapes
        meanShape = meanShape + alignedShapes{i};
    end
    meanShape = meanShape / numShapes;

    % Compute the Procrustes distance between the old and new mean shape
    [d, ~] = my_procrustes(reference, meanShape);
    fprintf('Left Iteration %d - Procrustes distance: %.6f\n', iter, d);

    % Check for convergence.
    if d < threshold
        fprintf('Left Convergence reached after %d iterations.\n', iter);
        reference = meanShape;
        break;
    end

    % Update the reference shape for the next iteration.
    reference = meanShape;
end

% Save final results for left hippocampus
final_mean_left = reference;
final_alignedShapes_left = alignedShapes; % final aligned shapes

% Compare Original and Mean Left Shape
redColors = repmat([1, 0, 0], size(initial_shape,1), 1);
blueColors = repmat([0, 0, 1], size(final_mean_left,1), 1); 
figure;
plot_mesh(initial_shape, face, struct('face_vertex_color', redColors));
hold on;
plot_mesh(final_mean_left, face, struct('face_vertex_color', blueColors));
title('Reference (red) and Mean (blue) Left Hippocampus Meshes');

% Overlay of All Aligned Left Shapes
figure;
hold on;
for i = 1:numShapes
    scatter3(final_alignedShapes_left{i}(:,1), final_alignedShapes_left{i}(:,2), final_alignedShapes_left{i}(:,3), 36, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
end
scatter3(final_mean_left(:,1), final_mean_left(:,2), final_mean_left(:,3), 60, 'r', 'filled');
title('Overlay of All Aligned Left Hippocampus Shapes with Final Mean (red)');
axis equal; 
grid on; 
view(3);

% Bar Plot of Procrustes Distances (Left)
distances_left = zeros(numShapes,1);
for i = 1:numShapes
    [d_i, ~] = my_procrustes(final_mean_left, final_alignedShapes_left{i});
    distances_left(i) = d_i;
end
figure;
bar(distances_left);
title('Procrustes Distance from Final Mean Shape for Each Left Hippocampus');
xlabel('Shape Index'); 
ylabel('Procrustes Distance');
grid on;

% PCA of Aligned Left Hippocampus Shapes
% Flatten each shape into a vector.
X_left = zeros(numShapes, numel(final_mean_left));
for i = 1:numShapes
    X_left(i,:) = final_alignedShapes_left{i}(:)';
end
[coeff_left, score_left, latent_left] = pca(X_left);
figure;
scatter(score_left(:,1), score_left(:,2), 50, 'b', 'filled');
xlabel('PC 1'); ylabel('PC 2');
title('PCA of Aligned Left Hippocampus Shapes');
grid on;

% Calculate the correlation matrix Before Alignment
X_before = zeros(numShapes, numel(final_mean_left));
for i = 1:numShapes
    X_before(i,:) = shapes{i}(:)';  % Flatten shape i
end

corrMatrix_before = zeros(numShapes, numShapes);
for i = 1:numShapes
    for j = 1:numShapes
         r = corrcoef(X_before(i,:), X_before(j,:));
         corrMatrix_before(i,j) = r(1,2);
    end
end

% Calculate the correlation matrix After Alignment
X_after = zeros(numShapes, numel(final_alignedShapes_left{1}));
for i = 1:numShapes
    X_after(i,:) = final_alignedShapes_left{i}(:)';
end

corrMatrix_after = zeros(numShapes, numShapes);
for i = 1:numShapes
    for j = 1:numShapes
         r = corrcoef(X_after(i,:), X_after(j,:));
         corrMatrix_after(i,j) = r(1,2);
    end
end

% Display the correlation matrices with colors
figure;
subplot(1,2,1);
imagesc(corrMatrix_before);
clim([0 1]); 
colorbar;
title('Left Correlation Matrix Before Alignment');
xlabel('Hippocampus Index');
ylabel('Hippocampus Index');

subplot(1,2,2);
imagesc(corrMatrix_after);
clim([0 1]);  
colorbar;
title('Left Correlation Matrix After Alignment');
xlabel('Hippocampus Index');
ylabel('Hippocampus Index');

% Compute the histogram of the upper-diagonal elements
upperDiag_before = corrMatrix_before(triu(true(size(corrMatrix_before)), 1));
upperDiag_after  = corrMatrix_after(triu(true(size(corrMatrix_after)), 1));

figure;
edges = 0:0.1:1; 
hold on;
histogram(upperDiag_before, 'Normalization', 'probability', 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'BinEdges', edges);
histogram(upperDiag_after,  'Normalization', 'probability', 'FaceColor', 'red',  'FaceAlpha', 0.5, 'BinEdges', edges);
xlabel('Correlation Coefficient');
ylabel('Probability');
title('Histogram of Left Upper Diagonal Matrix Before vs. After Alignment');
legend('Before Alignment', 'After Alignment');
grid on;
hold off;


%% Process Right Hippocampus
% Select one shape to be the initial approximate mean shape
folder = 'path\to\folder';
files = dir(fullfile(folder, '*_right.vtk')); % Get files with _right.vtk
numShapes = numel(files);
shapes = cell(1, numShapes);

% Load each .vtk file and extract the vertex data (transposed to Nx3)
disp('Loading right hippocampus files...');
for i = 1:numShapes
    fileName = fullfile(folder, files(i).name);
    [vertex, face] = read_vtk(fileName, 0); % Read .vtk file
    shapes{i} = vertex'; % Transpose
end

% Set the initial approximate mean shape
initial_shape = shapes{1};
reference = initial_shape;

% Set convergence parameters for GPA
maxIter = 10; % Maximum number of iterations
iter = 0;
threshold = 1e-10;

% GPA Iteration for Right Hippocampus:
while iter < maxIter
    iter = iter + 1;
    alignedShapes = cell(1, numShapes);

    % Align each shape to the current reference using Procrustes
    for i = 1:numShapes
        [~, Z] = my_procrustes(reference, shapes{i});
        alignedShapes{i} = Z;
    end

    % Calculate the new approximate mean shape from the aligned shapes
    meanShape = zeros(size(reference));
    for i = 1:numShapes
        meanShape = meanShape + alignedShapes{i};
    end
    meanShape = meanShape / numShapes;

    % Compute the Procrustes distance between the old and new mean shape
    [d, ~] = my_procrustes(reference, meanShape);
    fprintf('Right Iteration %d - Procrustes distance: %.6f\n', iter, d);

    % Check for convergence.
    if d < threshold
        fprintf('Right Convergence reached after %d iterations.\n', iter);
        reference = meanShape;
        break;
    end

    % Update the reference shape for the next iteration.
    reference = meanShape;
end

% Save final results for right hippocampus
final_mean_right = reference;
final_alignedShapes_right = alignedShapes; % final aligned shapes

% Compare Original and Mean Right Shape
redColors = repmat([1, 0, 0], size(initial_shape,1), 1); 
blueColors = repmat([0, 0, 1], size(final_mean_right,1), 1);
figure;
plot_mesh(initial_shape, face, struct('face_vertex_color', redColors));
hold on;
plot_mesh(final_mean_right, face, struct('face_vertex_color', blueColors));
title('Reference (red) and Mean (blue) Right Hippocampus Meshes');

% Overlay of All Aligned Right Shapes
figure;
hold on;
for i = 1:numShapes
    scatter3(final_alignedShapes_right{i}(:,1), final_alignedShapes_right{i}(:,2), final_alignedShapes_right{i}(:,3), 36, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
end
scatter3(final_mean_right(:,1), final_mean_right(:,2), final_mean_right(:,3), 60, 'r', 'filled');
title('Overlay of All Aligned Right Hippocampus Shapes with Final Mean (red)');
axis equal; 
grid on; 
view(3);

% Bar Plot of Procrustes Distances (Right)
distances_right = zeros(numShapes,1);
for i = 1:numShapes
    [d_i, ~] = my_procrustes(final_mean_right, final_alignedShapes_right{i});
    distances_right(i) = d_i;
end
figure;
bar(distances_right);
title('Procrustes Distance from Final Mean Shape for Each Right Hippocampus');
xlabel('Shape Index'); 
ylabel('Procrustes Distance');
grid on;

% PCA of Aligned Right Hippocampus Shapes
% Flatten each shape into a vector.
X_right = zeros(numShapes, numel(final_mean_right));
for i = 1:numShapes
    X_right(i,:) = final_alignedShapes_right{i}(:)';
end
[coeff_right, score_right, latent_right] = pca(X_right);
figure;
scatter(score_right(:,1), score_right(:,2), 50, 'b', 'filled');
xlabel('PC 1'); ylabel('PC 2');
title('PCA of Aligned Right Hippocampus Shapes');
grid on;

% Calculate the correlation matrix Before Alignment
X_before = zeros(numShapes, numel(final_mean_right));
for i = 1:numShapes
    X_before(i,:) = shapes{i}(:)';  % Flatten shape i
end

corrMatrix_before = zeros(numShapes, numShapes);
for i = 1:numShapes
    for j = 1:numShapes
         r = corrcoef(X_before(i,:), X_before(j,:));
         corrMatrix_before(i,j) = r(1,2);
    end
end

% Calculate the correlation matrix After Alignment
X_after = zeros(numShapes, numel(final_alignedShapes_right{1}));
for i = 1:numShapes
    X_after(i,:) = final_alignedShapes_right{i}(:)';
end

corrMatrix_after = zeros(numShapes, numShapes);
for i = 1:numShapes
    for j = 1:numShapes
         r = corrcoef(X_after(i,:), X_after(j,:));
         corrMatrix_after(i,j) = r(1,2);
    end
end

% Display the correlation matrices with colors
figure;
subplot(1,2,1);
imagesc(corrMatrix_before);
clim([0 1]); 
colorbar;
title('Right Correlation Matrix (Before Alignment)');
xlabel('Hippocampus Index');
ylabel('Hippocampus Index');

subplot(1,2,2);
imagesc(corrMatrix_after);
clim([0 1]);  
colorbar;
title('Right Correlation Matrix (After Alignment)');
xlabel('Hippocampus Index');
ylabel('Hippocampus Index');

% Compute the histogram of the upper-diagonal elements
upperDiag_before = corrMatrix_before(triu(true(size(corrMatrix_before)), 1));
upperDiag_after  = corrMatrix_after(triu(true(size(corrMatrix_after)), 1));

figure;
edges = 0:0.1:1; 
hold on;
histogram(upperDiag_before, 'Normalization', 'probability', 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'BinEdges', edges);
histogram(upperDiag_after,  'Normalization', 'probability', 'FaceColor', 'red',  'FaceAlpha', 0.5, 'BinEdges', edges);
xlabel('Correlation Coefficient');
ylabel('Probability');
title('Histogram of Right Upper Diagonal Matrix Before vs. After Alignment');
legend('Before Alignment', 'After Alignment');
grid on;
hold off;
