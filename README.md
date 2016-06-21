# Simulations to investigate the fused lasso for "Approximate Recovery in Changepoint Problems, from \ell_2 Estimation Error Rates"
The code to run the simulations and plot all the figures in the named paper.

# Dependencies
The follow R packages are required (listed in "source_header.R")
 - genlasso: For computing the fused lasso estimate
 - assertthat: For assertations to maintain code correctness
 - doMC: For parallel computation
 - foreach: For parallel computation

# How to run:
The scripts to run are in the "main.simulation" folder with the prefix "final". 
There are three seperate simulation pipelines.
 - Computing the fused lasso for a wide range of n values: Run "final_simulation.R" 
 (to generate a file called "final...RData") to be used in "fig3_filter-example.R",
  "fig4_haus.R" and "fig7_rates.R". 
 - Computing the fused lasso for n=774 and exploring the reduced filter: 
 Run "final_simulation_roc.R" (to generate a file called "final-ROC...RData") followed
 by "final_simulation_roc2.R" (to generate a file called "final-ROC2...RData")
 to be used in "fig5_filter-dist.R" and "fig7_rates.R"
 - Computing the fused lasso over the 2d graph:
 Run "final_2dgraph.R" (to generate a file prefixed "graph-dim") followed
 by "final_2dgraph_part2.R" (to generate a file called "graph-dim-...reduced.RData")
 to be used in "fig8_graph.R".

The remaining figures, generated in "fig1_intersim.R" and "fig2_data.R", can be
generated independently.

