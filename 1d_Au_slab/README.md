
    <script type="text/javascript"
            src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>

# One dimensional example
Calculation of the transmittance through a slab of gold modelled as a Drude-critical point medium. A comparison between FDTD code and the analytical solution is given.  

A new dispersive model for metals has been recently proposed that includes a single Drude term and $N$ critical point pairs . It has been proved in several works that the Drude-critical point (DCP) model can more accurately describe the dielectric dispersion of metals in infrared and optical frequencies than the widely used Drude-Lorentz medium . The relative dielectric permittivity of the DCP model, assuming $e^{j \omega t}$ time dependence, is described by
\begin{equation}
\begin{split}
\varepsilon (\omega) & = \varepsilon_{\infty} + \frac{\omega_D^2}{\omega(j
\gamma - \omega)} + \\ & \sum_{p=1}^{N} A_p \Omega_p \left( \frac{e^{j \phi_p}}{\Omega_p +
\omega - j \Gamma_p} + \frac{e^{- j \phi_p}}{\Omega_p - \omega + j
\Gamma_p} \right)
\end{split}
\end{equation}
where $\varepsilon_{\infty}$ is the relative permittivity at infinite frequency, $\omega_D$ and $\gamma$ are the plasma frequency and damping coefficient, respectively, of the Drude model, whereas $A_p$ is the amplitude, $\phi_p$ the phase, $\Omega_p$ the frequency and $\Gamma_p$ the broadening parameter, respectively, of the $p$-th critical point oscillator.

It can be observed that both the Drude and the CP contributions can be included in a more general term, which is an explicit function of the variable $j \omega$
\begin{equation}
\varepsilon (\omega) = \varepsilon_{\infty} +  \sum_{p=1}^{M} \frac{a_{1p} j \omega + a_{0p}}{b_{2p}(j \omega)^2 + b_{1p} j \omega + b_{0p}}
\label{er_general}
\end{equation}
where $M=N+1$.

The FDTD code is based on the following papers:

***K. P. Prokopidis, D. C. Zografopoulos and E. E. Kriezis***,  “Rigorous broadband investigation of liquid-crystal plasmonic structures using finite-difference time-domain dispersive-anisotropic models,” J. Opt. Soc. Am. B,  Vol. 30, No. 10, pp 2722-2730,  October 2013. [DOI:10.1364/JOSAB.30.002722](http://dx.doi.org/10.1364/JOSAB.30.002722)


***K. P. Prokopidis and D. C. Zografopoulos***, “A Unified FDTD/PML Scheme Based on Critical Points for Accurate Studies of Plasmonic Structures,” Journal of Lightwave Technology, vol. 31, pp. 2467-2476, Aug. 2013. [DOI: 10.1109/JLT.2013.2265166](http://dx.doi.org/10.1109/JLT.2013.2265166)
