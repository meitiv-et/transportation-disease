# The role of transportation in the spread of infectious disease

This repository contains C++ code to simulate a complex agent-based
model and visualize the results.  For the complete model description
and the summary of the findings please see [CARTEEH DataHub Data Story](https://carteehdata.org/library/document/the-role-of-transportatio-d7a7)

### Running the model:
1. Install `g++`, `gnuplot`, `ffmpeg`, and `make`
```bash
$ sudo apt-get install g++ gnuplot ffmpeg make
```
2. Compile the c++ code
```bash
$ make
```
3. Run the code (clean up the existing frames from a previous run first)
```bash
$ rm -f frame_*
$ ./two_component_grid_death numPeopleX numPeopleY numVehicles tau1
tau2 tau3 personContactRate vehicleDemandRate disinfectionRate
deathRate lambda
```
The meaning of the command line parameters are summarized in the table
below:

| Variable | Description |
| --- | --- |
| numPeopleX | people grid size in the X direction |
| numPeopleY | people grid size in the Y direction |
| numVehicles | number of vehicles |
| tau1 | number of days from infection onset to peak infectivity |
| tau2 | number of days from infection onset to recovery |
| tau3 | half-life of vehicle infectivity (in days) |
| personContactRate | healthy person to person contacts per day per person |
| vehicleDemandRate | healthy person to vehicle contacts per day per person |
| disinfectionRate | number of disinfections per day per vehicle (could be < 1) |
| deathRate | probability of death per day at the infection peak |
| lambda | parameter that determines transmission probability (see the model description of details) |
4. By default the executable generates the frames depicting the
   infectivity state of every person and vehicle `frame_0000.png`,
   etc.  A more detailed output can be enabled by un-commenting on or
   both `system.output` statements in the source code.
5. Render the frames using `gnuplot`
```bash
$ for f in frame_*gp; do gnuplot $f; done
```
6. Assemble the rendered frames into a movie using `ffmpeg`
```bash
$ ffmpeg -r 30 -f image2 -start_number 1 -i frame_%04d.png -vframes 1000 -vcodec libx264 -crf 25 -pix_fmt yuv420p epidemic.mp4
```

