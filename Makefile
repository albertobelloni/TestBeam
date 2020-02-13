SHELL := /bin/sh # sh is the default
.SILENT: clean slimprep plotclean

PHONY: clean

clean:
	mkdir -p results
	if ls roofit_e* >& /dev/null; then mv roofit_e* results/; \
	else echo "No roofit files"; fi
	if ls canvas_e* >& /dev/null; then mv canvas_e* results/; \
	else echo "No canvas files"; fi
	if ls canvas_gaussfit* >& /dev/null; \
	then mv canvas_gaussfit* results/; \
	else echo "No canvas Gauss fit files"; fi

libclean:
	rm -f *_C.d *_C.so *_C_ACLiC_dict_rdict.pcm

slimprep:
	mkdir -p Alignment_Plots Crud_test Denominator_Plots Energy_Plots \
	Noise_Plots \
	Original_Images/Efficiency_Maps_2D \
	Original_Images/Efficiency_Maps_X/No_Crud/Special_Bins \
	Original_Images/Efficiency_Maps_Y/No_Crud/Special_Bins \
	Overlayed_Plots Pedestal_Plots \
	Rotated_Images/Efficiency_Maps_2D/CMB_Plots \
	Rotated_Images/Efficiency_Maps_X/No_Crud/Special_Bins \
	Rotated_Images/Efficiency_Maps_Y/No_Crud/Special_Bins \
	Time_Slice_Plots
	echo "Created directory structure to hold plots!"


plotclean:
	rm -rf Alignment_Plots Crud_test Denominator_Plots Energy_Plots \
	Noise_Plots Original_Images Overlayed_Plots Pedestal_Plots \
	Rotated_Images Time_Slice_Plots
	echo "Deleted directory structure with sliman2017 plots"
