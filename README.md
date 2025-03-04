# Procrustes-Based Hippocampus Shape Alignment

This repository contains MATLAB scripts for hippocampus shape alignment using **Procrustes analysis**. The pipeline standardizes hippocampal meshes by removing translational, rotational, and scale differences, allowing for precise comparisons across subjects. This methodology is particularly relevant for neurodegenerative disease research, such as detecting structural changes in Alzheimer's disease.

## 📁 Project Structure

├── find_cc.m # Identifies and separates left and right hippocampus ├── read_vtk.m # Loads VTK mesh files ├── write_vtk.m # Saves processed hippocampus meshes in VTK format ├── check_face_vertex.m # Ensures vertex and face matrices have the correct format ├── my_procrustes.m # Custom pairwise Procrustes alignment ├── two_pair_registration.m # Registers left and right hippocampus shapes ├── hippocampus_GPA.m # Generalized Procrustes Analysis (GPA) ├── plot_mesh.m # Visualizes hippocampus meshes ├── step1_batch_separate_hippocampus.m # Automates hippocampus segmentation ├── step2_example_procrustes.m # Example script for pairwise Procrustes ├── step3_general_procrustes.m # Main script for GPA-based shape alignment └── README.md # Project documentation


## 🛠 Installation & Dependencies
This project is implemented in **MATLAB** and requires the following toolboxes:
- **Parallel Computing Toolbox** (for performance optimization)
- **Image Processing Toolbox** (for visualization support)

