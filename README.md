========================================================================
   PW_EIB_RING_BEAM_PARALLEL: Electron-Ion-Beam (E-I-B) Plasma Dispersion Solver
========================================================================

[Authors]      Yasuhito Narita (1) and Uwe Motschmann (1)

[Contact]     y.narita@tu-braunschweig.de

[Affiliation] (1) Institute of Theoretical Physics, Technical University of Braunschweig, Braunschweig, Germany.

[License]     MIT License

[DOI]         TBD

------------------------------------------------------------------------
1. OVERVIEW
------------------------------------------------------------------------
This repository hosts a self-contained scientific script for analyzing plasma instabilities.
The code **`pw_eib_ring_beam_parallel.py`** solves the complex dispersion relations for right-hand modes propagating parallel to the background magnetic field in an electron-ion-beam (e-i-b) plasma system. The word "pw" stands for plasma waves. The code operates within the displacement-current-free approximation, enforces charge and current neutrality conditions, and tracks growth rates ($\gamma$) for unstable modes induced by a ring velocity component.

The solver outputs production-ready publication figures in PDF format and serializes raw multidimensional calculation matrices into NumPy binary archive tracks (`.npz`) for downstream post-processing and analysis.

Details of the algorithm are presented in the following journal article:

 Narita, Y., and Motschmann, U.:
    Non-resonant beam instability driven by the ring-beam ion distribution,
    AIP Advances 16, in press (2026)

------------------------------------------------------------------------
2. REQUIREMENTS
------------------------------------------------------------------------
- Python 3.x
- NumPy
- Matplotlib (optional, for plotting)

------------------------------------------------------------------------
3. INSTALLATION
------------------------------------------------------------------------
The code pw_eib_ring_beam_parallel.py can run on the terminal.

------------------------------------------------------------------------
4. BASIC USAGE
------------------------------------------------------------------------
Set the following parameters: 

 (1) Particle mass parameters.
 The electron mass is set to unity. The bulk ion mass
 is set by the parameter omega_ce (with the minus sign),
 and the beam ion mass by omega_cb. 
 The default is the proton mass for both the bulk ions and the beam ions,
  omega_ce = -1836, omega_ci = 1, omega_cb = 1.

 (2) Particle density parameters.
 The bulk ion density (number density) is set by 
  the ion plasma frequency omega_pi (default 1000),
  and the beam density by omega_pb (default 100).
  The plasma frequency is proportional to the square density
  of the respective species. The electron density omega_pe is
  automatically determined by the charge neutrality condition.

 (3) Particle velocity parameters.
  Electron bulk speed v_e parallel to the magnetic field
  is set to zero (electron rest frame).
  Beam speed has two parameters. The parallel beam velocity
  v_b is scaled to the Alfven speed (default 20).
  The perpendicular beam velocity v_perp (ring distribution in the
  velocity space) is scaled to the Alfven speed, to (default 10).
  The bulk ion velocity is automatically computed from
  the current-free condition.

 (4) Wavenumber array (parallel wavenumbers).
  Give the starting wavenumber (can be negative values) and
  the ending wavenumber in units of the ion inertial wavenumer, 
  and the number of wavenumber arrays. Default is 
  kcw01_min = 0.2, kcw01_max = 0.2, and kcw01_ticks = 4000.


Run the code on the terminal by typing, e.g.,

\>\> python3 pw_eib_ring_beam_parallel.py

The outputs are the PDF file showing the dispersion relation
and the growth rate, and the NPZ binary data file containing
the wavenumbers and frequencies.


------------------------------------------------------------------------
5. CITATION
------------------------------------------------------------------------
If you use this code in your research, please cite it as:

 Narita, Y., and Motschmann, U.:
    Non-resonant beam instability driven by the ring-beam ion distribution,
    AIP Advances 16, in press (2026)

------------------------------------------------------------------------
6. LICENSE
------------------------------------------------------------------------
This project is licensed under the MIT License.
