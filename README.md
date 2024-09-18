# SecondaryStructureProtein
This code processes a molecular dynamics trajectory file to compute and visualize the secondary structure of a protein over time. It uses the mdtraj library to load a trajectory and restricts the analysis to protein atoms. The secondary structure of the protein is computed using the DSSP method, and the results are converted into integers representing different structural elements (e.g., alpha-helix, beta-sheet). A color map is created to visually distinguish these structures, and a plot is generated using matplotlib to display the secondary structure of the protein across time and residues. The plot is saved as "dssp_chains.png".
