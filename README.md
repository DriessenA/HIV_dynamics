# HIV_dynamics
Investigating slow CD4+ T cell depletion during chronic HIV infection. ODE model of memory CD4+ T cell population during HIV infection in which infection occurs upon cellular division through an infection chance. We also investigate the role of immune activation and homeostatic activation.

The Loop.R script contains the model and most analyses described in the manuscript (the memory CD4+ T cells effectively create a loop). The comp_HA.R script is to compare multiple function for scaling homeostatic activation (linear, convex, concave). 

This project was performed at Utrecht University under supervision of Prof. Rob de Boer, department Theoretical Biology and Bioinformatics. This project was an internship of the master Bioinformatics and Systems Biology at the Vrije Universitieit  and University of Amsterdam (supervised by examiner dr. ir. Anton Feenstra).

## Dependencies
Grind.R is needed for ODE simulations: http://tbb.bio.uu.nl/rdb/grindR/
