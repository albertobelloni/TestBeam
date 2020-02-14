# TestBeam
ROOT macros to analyze H2 testbeam ntuples

## Produce (most of) the plots and create the file with the energy tree

The main macro is +sliman2017.C+. First, start creating the directory structure
to hold the histograms (.C, .root, .pdf, and .png):

```
make slimprep
```

Then, run the main macro:

```
[] .L sliman2017.C++
[] doMaps()
[] doAlignmentPlots()
[] doEnergy()
[] doTimeSlices()
```

The file +energy_hists.root+ shall also be created. It contains a tree with
the energy of fiducial hits for each of the scintillator tiles.

## Run the fitter(s)

### Physics-driven

The physics-driven fitter requires the installation of Yi-Mu's RooFit PDFs.
Instructions can be found in the top of +fitter2017.C+
Once compiled, run the fits on condor:


```
condor_submit condorfits.jdl
```

The output is a set of roofit_?_?_?.root files. The first ? indicates whether
the fit used a binned histogram or a tree as input; the second ? indicates the
tile; the third ? is a set of flags (baseline/afterpulsing/full fit;
binned/unbinned fit; rebin value - used only if binned fit).

One shall then produce histograms, from the roofit_* files. This is
accomplished by running the +make_fitter_plots.C+ macro. Again, instructions
are provided in the top part of the macro.

Together with the plots, a text file with some fit results, +text_results.txt+,
is saved. Let us now run the multi-Gaussian fitter, and collect in the same
file the results of the other two <p.e.> estimators, before copying it to
the DN-18-007 directory.

### Multi-Gaussian

The main macro is +peakfinder.C+. Instructions on how to run it are reported
in the top part of that file. Once run, the +text_results.txt+ file will
also contain the results of the <p.e.> estimation using the histogram integral
and multi-Gaussian fit methods. It is now ready to be copied to the DN-18-007
directory.

## Collect the plots for the paper and the presentation

This is a bit weird, and I shall update the name of the target. As of now,
at this point we do:

```
make clean
make packfigs
```

The first command will move all files produced by the physics-driven fit to
a directory called +results+. The second command will scour all directories
and collect in a tarball, +dn-18-007_figs.tar+, the files needed to the paper
and the presentation.

## Clean all

The following command will remove +text_results.txt+, +dn-18-007_figs.tar+,
all the compiled code (_C.so, _C.d, and dictionary files), and the directories
created by the +make slimprep+ command. Directories containing the results from
the condor jobs are preserved.