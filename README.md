[CRes (Crater Resonance)](https://github.com/leighton-watson/CRes) is a one-dimensional (1D) numerical method for solving the linear acoustic wave equation within a volcanic crater. For a specified crater geometry and excitation source at the base of the crater, CRes computes the velocity and pressure at the crater outlet and can propagate the signal to an infrasound receiver some distance from the outlet. The linear acoustic wave equation is written in terms of two first-order differential equations for pressure and acoustic flow and is solved by CRes using a finite-difference frequency-domain method. CRes is written in MATLAB and runs efficiently on a standard desktop/laptop computer.

This repository contains data files and custom scripts relating to several papers that have used CRes for volcano crater resonance:
* Watson, L. M., Dunham, E. M., and Johnson, J. B. (2019) Simulation and inversion of harmonic infrasound from open-vent volcanoes using an efficient quasi-1D crater model, Journal of Volcanology and Geothermal Research, [https://doi.org:10.1016/j.jvolgeores.2019.05.007](https://doi.org/10.1016/j.jvolgeores.2019.05.007).
* Watson, L. M., Johnson, J. B., Sciotto, M., and Cannata, A. (2020) Changes in crater geometry revealed by inversion of harmonic infrasound observations: 24 December 2018 eruption of Mount Etna, Italy, Geophysical Research Letters, [https://doi.org/10.1029/2020GL088077](https://doi.org/10.1029/2020GL088077)
* Sciotto, M., Watson, L. M., Cannata, A., Cantarero, M., De Beni, E., Johnson, J. B. (2022) Infrasonic gliding reflects a rising magma column at Mount Etna (Italy), Scientific Reports, [https://doi.org/10.1038/s41598-022-20258-9](https://doi.org/10.1038/s41598-022-20258-9)
* Watson, L. M., Amato, G., Cannata, A., De Beni, E., Proietti, C., Sciotto, M. (submitted) Long-Lasting Infrasonic Gliding with Multiple Harmonics at Mount Etna, Volcanica.

  The main modeling code is found at [https://github.com/leighton-watson/CRes](https://github.com/leighton-watson/CRes). Make sure to install this and point to that directory before trying to run the scripts from this repository.

  ### Who do I talk to? ###
* Leighton Watson: leighton.watson@canterbury.ac.nz
