# Background
The use of renewable energy sources, such as solar, wind and hydraulic energies, is very old,
they have been used since many centuries before our time, by example, mills use wind energy to
grind wheat. Its applications ceased during the "Industrial Revolution", at which time, due to the
low price of petroleum, they were abandoned. During recent years, due to the increase in fossil
fuel prices and the environmental problems caused by the use of conventional fuels, we are
reverting back to renewable energy sources.

Renewable energies are virtually inexhaustible, clean and they can be used in a decentralized
way. One mayor problem of this kind of energy is the randomness and variability of its
production because of its dependence on meteorological variables, which are hard to predict.
Among them, wave energy has the advantage of being relatively concentrated, constant and
predictable. Unlike other renewal energies such as solar energy, which can only by used during
the daylight, wave energy is usable throughout the whole day, but it‟s only available on the sea.
The study and prediction of the wave energy available in a certain area is essential for its proper
viability, performance and design of systems able to transform that energy into usable energy.

In this project, Soft-Computing techniques are proposed for the prediction of wave
characteristics in a sea area, using sea wave parameters from other places, in the same period of
time. Neural networks are used to predict the height of the sea waves using other parameters as
input; and metaheuristic search techniques are used to find the subset of these parameters that
produces the best prediction.

# Implementation

The project is completely implemented in C using the CBLAS library and POSIX threads for optimizing the speed.
The main problem is an Feature Selection problem where we have to select the variables that can better predict
the height of a wave from a set of variables. This is done by combining a fast regression algorithm, the 
Extreme Learning Machine and Discrete Optimization Local Search Heuristics sucha as Simmulated Anneling 
and Evolutionary Algorithms.

Codes are comments are fully in English but the documentation to use the tool and the thesis describing the
thesis are in Spanish since this was my final thesis of my Telecommincations Engineering Degree. There is fast and nice
implementation of the Penrourse-Moore inverse using C which can be reused, as well as some interesting algorithms.

The code has been tested with the valgrind tool to ensure that no memmory is left without allocation.
In order to build the program follow the Makefile, everything is explained in the Thesis that can be found inside the "Documentation" folder.
