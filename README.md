# Analysis of the spatial ecology of siphonophores in the Offshore Central California Ecosystem

This pipeline uses real positions in Z (depth), relative abundances, and patchiness (degree of spatial aggregation, when available from transects) from the VARS data, and generates random X and Y components over thousands of iterations. In each iteration, my pipeline will record the relative frequencies of 2 point objects being within threshold (i.e. 1m) distance for each siphonophore-prey taxon couple. The output is a relative encounter rate matrix for each siphonophore species (rows) and prey taxon (columns) in a scale from 0 to 1. Adjustments on the threshold distance can be made for different siphonophore and prey species, accounting for siphonophore colony length, tentacle length, prey motility, and prey size.

## Goal:

Simulate and compare 3D spatial distribution patterns of different siphonophore species and other zooplankton (potential prey) from ROV annotations' bathymetric summary statistics and patchiness metrics.

Calculate distance distributions between different planktonic taxa.

Calculate relative encounter probabilities for each pair of siphonophore-prey species.

Estimate an expected dietary covariance matrix for siphonophore species pairs given their distributions relative to prey.

## Spatial Analysis in 2D

Generate PPP objects using sp package.
Calculate spatial sumary statistics: Spherical contact

Take a 500 points sample of the large XYZ matrix, calculate euclidean pairwise distances between species (from the set of points that belong to each species using a unique identifier loop), plot heatmap.
