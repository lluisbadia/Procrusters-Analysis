# Procrustes-Based Hippocampus Shape Alignment

This repository contains MATLAB scripts for hippocampus shape alignment using **Procrustes analysis**. The pipeline standardizes hippocampal meshes by removing translational, rotational, and scale differences, allowing for precise comparisons across subjects. This methodology is particularly relevant for neurodegenerative disease research, such as detecting structural changes in Alzheimer's disease.

## 📁 Project Structure

find_cc.m → Identifies and separates the left and right hippocampus from a given mesh.
read_vtk.m → Loads VTK mesh files containing hippocampus segmentation data.
write_vtk.m → Saves processed hippocampus meshes in VTK format.
check_face_vertex.m → Ensures the vertex and face matrices have the correct dimensions.
my_procrustes.m → Custom implementation of pairwise Procrustes alignment.
two_pair_registration.m → Registers left and right hippocampus shapes across subjects.
hippocampus_GPA.m → Performs Generalized Procrustes Analysis (GPA) to compute an average hippocampus shape.
plot_mesh.m → Visualizes hippocampus meshes before and after alignment.
step1_batch_separate_hippocampus.m → Automates hippocampus segmentation for multiple files.
step2_example_procrustes.m → Example script for pairwise Procrustes registration.
step3_general_procrustes.m → Main script for GPA-based shape alignment.


## 🛠 Installation & Dependencies
This project is implemented in **MATLAB** and requires the following toolboxes:
- **Parallel Computing Toolbox** (for performance optimization)
- **Image Processing Toolbox** (for visualization support)

