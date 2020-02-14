# TestBeam

ROOT macros to analyze H2 testbeam ntuples

## Produce (most of) the plots and create the file with the energy tree

The main macro is _sliman2017.C_. First, start creating the directory structure
to hold the histograms (.C, .root, .pdf, and .png):

```
make slimprep
```

Note the location of the H2 test beam ntuples. They currently are in:

```
LXPLUS: /afs/cern.ch/work/a/abelloni/CERN_July2017_TestBeam/SLIM_2/
HEPCMS: /data/users/abelloni/CERN_TB_Jul17/SLIM_2/
```

One can explicitly indicate the file directory when using the commands below,
or simply uncomment the relevant line in _sliman2017.h_ to conveniently
select the proper default value for the parameter `slim_dir`.
Then, run the main macro:

```
[] .L sliman2017.C++
[] doMaps()
[] doAlignmentPlots()
[] doEnergy()
[] doTimeSlices()
```

The file _energy_hists.root_ shall also be created. It contains a tree with
the energy of fiducial hits for each of the scintillator tiles, and binned
energy distributions. This file is used as input to the fitters.

## Run the fitter(s)

### Physics-driven

The physics-driven fitter requires the installation of Yi-Mu's RooFit PDFs.
Instructions can be found in the top of _fitter2017.C_.
Once compiled, run the fits on condor:


```
condor_submit condorfits.jdl
```

The output is a set of _roofit\_?\_?\_?.root_ files.
The first ? indicates whether the fit used a binned histogram or a tree
as input; the second ? indicates the tile;
the third ? is a set of flags (baseline/afterpulsing/full fit;
binned/unbinned fit; rebin value - used only if binned fit).

One shall then produce histograms, from the _roofit\_\*_ files. This is
accomplished by running the _make\_fitter\_plots.C_ macro. Again, instructions
are provided in the top part of the macro.

Together with the plots, a text file with some fit results, _text\_results.txt_,
is saved. Let us now run the multi-Gaussian fitter, and collect in the same
file the results of the other two <p.e.> estimators, before copying it to
the DN-18-007 directory.

### Multi-Gaussian

The main macro is _peakfinder.C_. Instructions on how to run it are reported
in the top part of that file. Once run, the _text\_results.txt_ file will
also contain the results of the <p.e.> estimation using the histogram integral
and multi-Gaussian fit methods. It is now ready to be copied to the DN-18-007
directory.

## Collect the plots for the paper and the presentation

At this point we do:

```
make fitclean
make packfigs
```

The first command will move all files produced by the physics-driven fit to
a directory called _results_. To be honest, I think it is not necessary, since
the second command will scour all directories and collect in a tarball,
_dn-18-007\_figs.tar_, the files needed to the paper and the presentation.
The directory does look nicer, if one runs `make fitclean`.

## Clean all

The following command will remove _text\_results.txt_, _dn-18-007\_figs.tar_,
all the compiled code (\_C.so, \_C.d, and dictionary files), and the directories
created by the `make slimprep` command:

```
make cleanall
```

Directories containing the results from the condor jobs are preserved.