# 2015_GECCO_Bongard
Contains all of the material for replicating the paper entitled 'Evolving robot morphology
facilitates the evolution of neural modularity and evolvability'

The paper itself is available here: http://www.cs.uvm.edu/~jbongard/papers/2015_GECCO_Bongard.pdf

This readme file will help you to replicate the results from the paper. This is done in a gradual
fashion: each step enables you to replicate one of the figures or tables in the paper.

Step 1: Replicate Figure 1.

- Click on Fig1.pdf in this directory.

Step 2: Replicate Figure 2.

- This figure shows an evolved robot from experiment set 3 in Table 1. This robot was produced
  by allowing both the neural network and morphology of the robot to evolve. A bi-objective
  fitness function is used that only selects for age and grasping (the AFPO optimization
  method was used here).

- Enter the Evolve_For_Age_And_Grasping directory.

- Compile the C++ code by running ./makeModularity

- This produces an executable file called Modularity.

- Run this code at the command line with the following three arguments:

- Modularity 0 1 1

- This tells the code to run using random seed 0 (the first argument);

- to select for robot controllers that settle into a fixed attractor (the second argument); and

- to evolve the robot's morphology along with its controller (the third argument).

- As the file runs, it will print out the percent of progress that has been achieved so far.
  When progress reaches 100%, the program will finish.

- If the code is taking too long on your computer, open constants.h and reduce MAX_GENERATIONS.
  Recompile and re-run.

- When the code finishes, results are stored in the Data/ directory.

- Open Data/results_0.txt to see the results. Each row reports the modularity (m) and grasping
  ability (g) of the best robot in the population up to that point. Disregard the other data
  for now.

- Now let's visualize the behavior of the best robot found by the end of your run. To do so,
  copy the Robot_Matrix files found in Data/ into the Visualization/ directory.

- We will use Python, NumPy, SciPy and MatplotLib to do the visualization. If you are missing any
  of these, install them now.

- In Visualization/ type 'python Robot_Draw.py'. This will re-create something similar to Figure 2.

- Your visualization may not contain all four panels. This indicates that, by the end of your run,
  no robot had succeeded in all four environments. If this is the case, go back to constants.h
  and extend the length of the run by increasing MAX_GENERATIONS, re-compiling, and re-running. 
