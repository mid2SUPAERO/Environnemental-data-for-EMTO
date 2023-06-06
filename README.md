# Environnemental data for EMTO (EEMTO)

The Environnemental data for EMTO (EEMTO) aims to find the most ecological way to produce a work-piece which is to say that producing impact on environment is the lowest between proposed processes.

The base of EEMTO was Efficiency Multiscale Topology Optimisation (EMTO) by E.Duriez.  You can find initial code here: https://github.com/mid2SUPAERO/EMTO.

EEMTO is also based on the code Code Pando MultiStart that affords to find the lowest compliance based on different initial designs of structure. However, this code also affords to chose just one initiale design. In the majourity of calculation that we did we used just one initiale design with a view to not take up a lot of computer memory.

# Data base EEMTO 

EEMTO needs to have a material database that has been created using soft "Ansys EduPACK", extensions "Level 3". This database is created in Microsoft Excel. This is a good platform to work in differents programms (e.g. MatLab or others). 

In this data base there are information about CO2 footprint (production, recycling), heat energy (for polymers), energy consumtion for couple "material+fabrication process" Young's modulus, Poisson's ratio, density.

Furthermore there are Young's modulus ratio between conventional proceeds and additive manufacturing proceeds. Indeed, Young's modulus value depends on proceed method and it's crucially important to know the real mechanical property value. This permits to achive correct result (to find a material with good balance between mechanical properties and proceed energy consumption).

Although we didn't consider special requirements for our MBB case that can disbalance our choice "finding a material with good balance between mechanical properties and proceed energy consumption". These special requirements are resistance to corrosion (naval application), resistance to high temperatures (high speed aircrafts, recoverable spacecraft) etc. Due to their special requirements they need to have special mechanical propreties and the ecological question is not so much important. 

# References

[1] Duriez, E., Morlier, J., Charlotte, M., & Azzaro-Pantel, C. (2021). A well connected, locally-oriented and efficient multi-scale topology optimization (EMTO) strategy. Structural and Multidisciplinary Optimization, 1-24.

[2] Informations provenant du Logiciel Ansys EduPack. 2022

[3] Maëlle DUGUIN Léa MEIH-Sarah ZAYE Emilie BRUNO, Julie CAVAILLES. Impression 3d de materiaux legers et eco-conçus pour des applications en mobilite. 2021
