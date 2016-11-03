network SEIR model
============

Attempting to replicate the results from [Trapman et al (2016): Inferring R0 in emerging epidemics—the effect of common population structure is small](http://rsif.royalsocietypublishing.org/content/13/121/20160288)

Model
---------

* SEIR model
* Uses gillespie algorithm

Parameters
-----------

* Transmission rate ($\beta$): 0.02 days <sup>-1</sup>
* Mean duration of latent period ($\sigma$): 9.4 days
* Mean duration of infectious period ($\gamma$): 5 days
* Number of initially infected individuals ($I_0$): 10

Simulation
-----------

* Start with 10 initial infected individuals
    * Randomly chosen for each run