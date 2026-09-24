# musrfit-dynGssEALF
## WHAT IS THIS?

This set of source files implements the longitudinal-field muon spin relaxation function defined in [J. Phys. Soc. Jpn. 93, 044602 (2024)](https://journals.jps.jp/doi/10.7566/JPSJ.93.044602) as a user function for [musrfit](https://rmlmcfadden.github.io/musr/musrfit/).
This function describes muon spin relaxation in fluctuating local fields that exhibit the Edwards-Anderson-type autocorrelation for the nu_{muon} = 0 case in eq.11 of the above paper, parameterized with Delta, nu, Q, and B_{LF}.
In the current version, the function returns 0 if provided parameters are out of range (supported ranges: 0 <= Delta*t <= 6, nu >= 0, 0 <= Q <= 1, B_{LF} >= 0).
The source files for an experimental function for the nu_{ion}, nu_{muon} >= 0 case are also bundled in the dual_diffusion_model directory. See README in that directory for more information.

## AUTHOR

T. U. Ito  
Japan Atomic Energy Agency

## INSTALLATION

[Requirements] musrfit/ROOT running environment on any platform (see the [musrfit documentation](https://lmu.web.psi.ch/musrfit/user/html/index.html) for more information)

1. Unzip dynGssEALF_tbl_*.cpp.gz

2. Edit the "INSTBASE" variable in makefile to fit your installation location.
(By default, header and library files will be installed in $ROOTSYS/include/ and $ROOTSYS/lib/, respectively.)

3. Make and make install

## USAGE

The typical call through the msr-file would be
```
###############################################################
FITPARAMETER
#      Nr. Name        Value     Step      Pos_Error  Boundaries
        1 Alpha       1.         0.01        none        0        2     
        2 Asy         0.25       0.01        none        0.1     0.3
        3 LF          0.001      0.00        none        0        0.1
        4 Q           0.5        0.01        none        0.3     1
        5 Delta       0.2        0.01        none        0.15   0.25
        6 nu          1.5        0.01        none        0        10
##############################################################
THEORY
asymmetry      2
userFcn  libdynGssEALFLibrary.so   dynGssEALF   3   4   5   6 (LF Q Delta nu)
###############################################################
```
where the values for LF, Delta, and nu are given in Tesla, (microsec)^-1 and MHz, respectively.

## IMPORTANT NOTES ON USAGE

Please keep the following points in mind when using this function:

- **Primary Purpose:** This function was developed primarily for simulating $\mu$SR relaxation spectra. The fitting capability is implemented purely as an optional feature.
- **Reference & Intent:** For full details on the development background, physical intent, and examples of simulation applications, please refer to [arXiv:2505.09147](https://doi.org/10.48550/arXiv.2505.09147).
- **Parameter Correlation in Fitting:** When employing this function for parameter fitting, be aware that strong correlations exist among the parameters. It is strongly recommended to reduce the number of free parameters by applying appropriate physical constraints.
- **Disclaimer on Unconstrained Fitting:** Unconstrained fitting with all parameters left free provides no guarantee of yielding physically meaningful results. Users are responsible for evaluating and interpreting fit results within a physically sound framework.
- **Application to Ion Diffusion:** When applying this function to the analysis of ion diffusion, please exercise caution. As emphasized in Section *V. DISCUSSION AND SUMMARY* of [arXiv:2505.09147](https://doi.org/10.48550/arXiv.2505.09147), there are numerous essential physical factors that are not captured by this simplified model, requiring careful and critical interpretation of any analytical results.

## CITATION

T. U. Ito and R. Kadono, Distinguishing Ion Dynamics from Muon Diffusion in Muon Spin Relaxation,
[J. Phys. Soc. Jpn. 93, 044602 (2024)](https://journals.jps.jp/doi/10.7566/JPSJ.93.044602). DOI: 10.7566/JPSJ.93.044602


