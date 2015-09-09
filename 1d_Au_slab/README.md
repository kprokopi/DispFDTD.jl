# Transmission through film of gold
Various disperion models have been proposed such as the Drudeand the Drude-Lorentz with one or more poles to properly describe metals in the UV and optical range.

Recently, a new dispersive model has been proposed [1] that includes a single Drude term and N critical point pairs. It has been proved in several works that the Drude-critical point (DCP) model can more accurately describe the dielectric dispersion of metals in infrared and optical frequencies than the widely used Drude-Lorentz medium. 

Calculation of the transmittance through a 20-nm slab of gold modelled as a Drude-critical point medium is performed using the ADE-FDTD method proposed in [2,3].
In Julia REPL run 

    include("fdtd22_Drude_CP_slab_ADE.jl")
    include("fdtd22_air.jl")

Two files (total.txt and inc.txt) are saved on the disk. Then run the code
    
    include("Transmittance.jl")

to calculate the transmittance. A comparison between FDTD code and the analytical solution is calculated and is given in the following figure.

![Transmittance](Transmittance.png)  


#References
[1] ***A. Vial and T. Laroche***, “Description of dispersion of metals by means of the critical points model and applications to the study of resonant structures using FDTD method,” J. Phys. D:Appli. Phys., col. 40, pp. 7152-7158, 2007. 

[2] ***K. P. Prokopidis, D. C. Zografopoulos and E. E. Kriezis***,  “Rigorous broadband investigation of liquid-crystal plasmonic structures using finite-difference time-domain dispersive-anisotropic models,” J. Opt. Soc. Am. B,  Vol. 30, No. 10, pp 2722-2730,  October 2013. [DOI:10.1364/JOSAB.30.002722](http://dx.doi.org/10.1364/JOSAB.30.002722)


[3] ***K. P. Prokopidis and D. C. Zografopoulos***, “A Unified FDTD/PML Scheme Based on Critical Points for Accurate Studies of Plasmonic Structures,” Journal of Lightwave Technology, vol. 31, pp. 2467-2476, Aug. 2013. [DOI: 10.1109/JLT.2013.2265166](http://dx.doi.org/10.1109/JLT.2013.2265166)
