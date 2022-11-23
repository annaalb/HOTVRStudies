Two small programs to create displays of jets.

make_displays creates plots in phi-eta space with colored areas to illustrate the behaviour of jet finding algorithms.
Create the executable with
```
make make_displays
```
The input and output files have to be provided in make_displays.C. An example input is given in files/example/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root.


make_history creates a large number of gifs, showing each step of the clustering in phi-eta space. In each step, two pseudo-jets are clustered or called a jet. This can be used to examine the behaviour of jet finding algorithms in detail.
Create the executable with
```
make make_history
```
The input and output files have to be provided in make_history.C. The code needs a file with the jet clustering history and a file with the input particles. Examples are given in files/example. In addition, the user can provide a file with jet display histograms from the make_displays program. These will be used to color the final jets in the same color to have matching pictures. An example is given in files/example/uhh2.AnalysisModuleRunner.MC.ttbar_pythia8_flat_nnpdf23.root

Have fun!
