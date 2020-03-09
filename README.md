# Final Project

This repository is dedicated to the final project of the discipline "Data Explorative Analysis".

The theme of the project is the analysis of the position of individuals in a model of Evolution and Speciation by Marcus Aguiar, and programmed by Irina Lerner, as a Scientific Initiation funded by FAPESP.

If you want to find the model code, you can access it at https://github.com/IriLrnr/EvIBM. There you can find the instructions for running the code.

The context for this work is a problem: during code execution, individuals are accumulating together in points of the space, but the previous implementation of the model predicts that they should be homogenicaly distributed. To solve the problem, it is necessary to study it and analyze the map of the individuals' positions, and other features of the program.

## The scripts

To visualize what is happenning with the location on individuals throughout the simulation, we need to compare different generations. To do this, I am creating animations, showing the change in position during runtime, using gganimate.

Note: To include library "gganimate", I had to manually install libxml-2.0.
