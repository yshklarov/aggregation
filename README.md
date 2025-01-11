Aggregation
===========

This is a toy statistical mechanics simulation of atoms randomly aggregating into large molecules. Initially, the grid
is populated with matter at random positions (each grid cell contains matter with probability 'density' as specified by
the user). If two adjacent cells contain matter, they are clumped together into a single cluster. Then, at every time
step, each cluster moves in one of four cardinal directions picked at random. As soon as two cluster touch, they merge
into one. The simulation ends once some pre-specified proportion (the 'threshold') of total matter belongs to a single
cluster.


Screenshot
----------

Here is a snapshot of a 1000x1000 simulation in progress.

![screenshot](screenshots/screenshot.png)


Results
-------

This plot shows the number of time steps required to achieve a large cluster of at least 90% density, as a function of
lattice size, with 10% density. Each data point is a mean of 10,000 trials.

![results](results/results_10000_trials_three_densities.png)


Running
-------

On Linux, run 'make' to build. To build on other platforms, refer to the
[fenster readme](https://github.com/zserge/fenster).
