# TestBeam

ROOT macros to analyze H2 testbeam ntuples

## Prepare slimmed ntuples

The macro `slimmer_2017.C` is used to produce ntuples with smaller
event footprint than the ntuples produced with whichever package will
produce the next test-beam ntuples. For the 2017 test-beam data, we
used https://github.com/BaylorCMS/HCALTestBeam (more details at the
bottom of this README.md file).

There is a small inconveniency here. We need to:

- make the slimmed ntuples once

- run the macro which produces the alignment plots, `sliman2017.C`
  (run the `doAlignmentPlots()` function)

- obtain the x- and y- offset of the wire chambers, as printed out at
  the end of step above

- write the offsets in the top part of `slimmer_2017.C`and
  `sliman2017.h`

- well, run again the slimming job...

## Produce (most of) the plots and create the file with the energy tree

The main macro is `sliman2017.C`. First, start creating the directory
structure to hold the histograms (`.C`, `.root`, `.pdf`, and `.png`):

```
make slimprep
```

Note the location of the H2 test beam ntuples. They currently are in:

```
LXPLUS: /afs/cern.ch/work/a/abelloni/CERN_July2017_TestBeam/SLIM_2/
HEPCMS: /data/users/abelloni/CERN_TB_Jul17/SLIM_2/
```

One can explicitly indicate the file directory when using the commands below,
or simply uncomment the relevant line in `sliman2017.h` to conveniently
select the proper default value for the parameter `slim_dir`.
Then, run the main macro:

```
[] .L sliman2017.C++
[] doMaps()
[] doAlignmentPlots()
[] doEnergy()
[] doTimeSlices()
```

The file `energy_hists.root` shall also be created. It contains a tree with
the energy of fiducial hits for each of the scintillator tiles, and binned
energy distributions. This file is used as input to the fitters.

## Run the fitter(s)

### Physics-driven

The physics-driven fitter requires the installation of Yi-Mu's RooFit PDFs.
Instructions can be found in the top of `fitter2017.C`.
Once compiled, run the fits on condor:


```
condor_submit condor_fits.jdl
```

Check the status of your condor jobs with:

```
condor_q
```

The output is a set of `roofit_?_?_?.root` files.
The first ? indicates whether the fit used a binned histogram or a tree
as input; the second ? indicates the tile;
the third ? is a set of flags (baseline/afterpulsing/full fit;
binned/unbinned fit; rebin value - used only if binned fit).

One shall then produce histograms, from the `roofit_*` files. This is
accomplished by running the `make_fitter_plots.C` macro. Again, instructions
are provided in the top part of the macro. An important note: the function
`make_all_fitter_plots()` will use the `roofit_*` files contained in the
`results` directory. They are copied there upon running `make fitclean`, but
it is easy enough to copy them by hand.

Together with the plots, a text file with some fit results,
`yield_results.txt`, is saved. Let us now run the multi-Gaussian fitter,
and collect in the same file the results of the other two <p.e.> estimators,
before copying it to the DN-18-007 directory.

### Multi-Gaussian

The main macro is `peakfinder.C`. Instructions on how to run it are reported
in the top part of that file. Once run, the `yield_results.txt` file will
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
a directory called `results`. To be honest, I think it is not necessary, since
the second command will scour all directories and collect in a tarball,
`dn-18-007_figs.tar`, the files needed to the paper and the presentation.
The directory does look nicer, if one runs `make fitclean`.

## Clean all

The following command will remove `yield_results.txt`, `dn-18-007_figs.tar`,
all the compiled code (`C.so`, `C.d`, and dictionary files), and the directories
created by the `make slimprep` command:

```
make cleanall
```

Directories containing the results from the condor jobs are preserved.

## Production of ntuples

This is a moderately complicate step. The code is available at
https://github.com/BaylorCMS/HCALTestBeam.

The branch `CMSSW_8_1_0_pre7` seems to be the one I used last.
Two notes about them:

- `UserCode/H2TestBeamAnalyzer/tb_ana.py`: three `ROOT.gROOT.ProcessLine`'s
  contain `UChar_t`'s, which I replaced with `Float_t`
- the EMap ultimately used may be `EMAP-14JUL2017_Phase1_RM1RM2_Phase2_RM3.txt`
  and `EMAP-22JUL2017_Phase2_RM3.txt`, which I put in this repository for the
  sake of completeness