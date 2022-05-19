# Introduction

Spinodal decomposition is the mechanism in which a single thermodynamic phase seprates into two without nucleation. This is often observed in binary alloy systems where a mixture of two metal species seperate into two coexisting phases. The time evoloution of a binary alloy system is commonly modelled using the Cahn-Hilliard equation, which is a non-linear parabolic partial differntial equation. This software solves the Cahn-Hilliard equation using either finite differnces and explicit time integration, or via spectral methods using fast fourier transforms. SpinDec2 also has the capability of alterning the definition of the chemical mobility of the system to be dependent on the tempearture or concentration fields or both. The definition for these mobility fields are described by the Darkens equation and Boltzmann statistics, and are outlined in Y. H. Wang et al., J. Appl. 12 85102 (2019)). SpinDec2 has also been parralleised using OpenMP and MPI to optimise the code to run with incresed efficincey and uses the netcdf package to write the output and checkpoint files.


# All input is given in the form 'key=value'
        
Concentration_max=&lt;real&gt;  ==&gt; upper bound of uniform distribution that initialises concentration field

Concentration_min=&lt;real&gt;  ==&gt; lower bound of uniform distribution that initialises concentration field

Domain_x_size=&lt;real&gt;  ==&gt; size of x space

Domain_y_size=&lt;real&gt;  ==&gt; size of y space

Mobility_A=&lt;real&gt;  ==&gt; Atomic mobility of species A (for temperature independent mobility)

Mobility_B=&lt;real&gt;  ==&gt; Atomic mobility of species B (for temperature independent mobility)

free_energy_gradient_parameter=&lt;real&gt;  ==&gt; strength of interfacial term (kappa)

Bulk_free_energy=&lt;real&gt;  ==&gt;  bulk free energy scale factor (A) for case f(c) = Ac^2(1-c)^2

Checkpointing_interval=&lt;intl&gt;  ==&gt; checkpointing interval

Checkpoint_output_file=&lt;str&gt;  ==&gt; checkpoint file name

f(c)=&lt;real, real, ...&gt;  ==&gt; list of polynomial coefficient values to set bulk free energy

Max_time=&lt;real&gt; ==&gt; End time of simulation

dF_tolerance=&lt;real&gt; ==&gt; Tolerence limit in change in total free energy

Random_Seed=&lt;real&gt; ==&gt; Seed to reproduce result

Use_input=&lt;real&gt; ==&gt;  overwrites checkpoint meta data when set to one

Excitation_A=&lt;real&gt; ==&gt;  Exciation energy of species A (for temperature dependent mobility)

Excitation_B=&lt;real&gt; ==&gt;  Exciation energy of species B (for temperature dependent mobility)

temperature_max=&lt;real&gt;  ==&gt; upper bound of uniform distribution that initialises temperature field

temperature_min=&lt;real&gt;  ==&gt; lower bound of uniform distribution that initialises temperature field

Problem=&lt;strl&gt;  ==&gt; defines which mobility definition to use ('Consatnt', 'NonTemp', 'Temp')

Stabilization_Term=&lt;int&gt;  ==&gt; controls the stability of spectral method solver


Authors: Anas Siddiqui, Ben Gosling, Dyaln Morgan, Geraldine Anis, Matyas Parrag
