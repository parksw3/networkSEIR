network SEIR model
============

Attempting to replicate the results from [Trapman et al (2016): Inferring R0 in emerging epidemics—the effect of common population structure is small](http://rsif.royalsocietypublishing.org/content/13/121/20160288)

Model
---------

* SEIR model
* Uses gillespie algorithm

Parameters
-----------

* Transmission rate (beta): 0.02 days <sup>-1</sup>
* Mean duration of latent period (sigma): 9.4 days
* Mean duration of infectious period (gamma): 5 days
* Number of initially infected individuals (I(0)): 10

Simulation
-----------

* Start with 10 initial infected individuals
    * Randomly chosen for each run
* Model is run until epidemic dies out
* Only doing 100 runs for now... (Takes about 8 hours)
    
Calculation
---------------

* Infection based generation
    * Seems sensitive on the definition of generation k
    * For now generation k is defined as the first generation that has at least 75 infected individuals
* Initial exponential growth
    * Fit a straight line (Current method)
    * Estimation based on two points?

Current result
----------

![figure 5](https://github.com/parksw3/networkSEIR/blob/master/figures/Trapman.png)