SHELL := /bin/sh # sh is the default
.SILENT: help fitclean libclean slimprep plotclean

.PHONY: help

help:
	echo ""
	echo "Here is what you can do:"
	echo ""
	echo " help:      print this message"
	echo ""
	echo " list:      print list of all possible targets"
	echo ""
	echo " fitclean:  move physics-driven and multi-Gaussian fit results"
	echo "            to results/; create results/ if it does not exist"
	echo ""
	echo " libclean:  remove all _C.so, _C.d, and dict files"
	echo ""
	echo " slimprep:  create directory structure to hold non-fit plots"
	echo ""
	echo " plotclean: remove the directory structure with non-fit plots"
	echo ""
	echo " packfigs:  scour all directories and save .png files for note"
	echo "            and presentation in a tarball"
	echo ""
	echo " cleanall:  clean directory to -almost- checkout state; does not"
	echo "            remove results/ directory"
	echo ""

list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | \
	awk -v RS= -F: '/^# File/,/^# Finished Make data base/ \
	{if ($$1 !~ "^[#.]") {print $$1}}' | grep -v PHONY | grep -v \.png | \
	sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

fitclean:
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
	echo "Removed all _C.d, _C.so, ACLic_dict_rdict.pcm files"

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

FIGS =  Overlayed_Finger_YEff_nbins.png \
	Overlayed_Finger_XEff_nbins.png \
	Overlayed_Sigmas_YEff_nbins.png \
	Overlayed_Sigmas_XEff_nbins.png \
	align_2.png \
	align_3.png \
	efficiency_map_EJ_260.png \
	efficiency_map_SCSN_81F2.png \
	efficiency_map_rotEJ_200.png \
	efficiency_map_rotEJ_260.png \
	efficiency_map_rotEJ_260_2P.png \
	efficiency_map_rotSCSN_81F1.png \
	efficiency_map_rotSCSN_81F2.png \
	efficiency_map_rotSCSN_81F3.png \
	efficiency_map_rotSCSN_81F4.png \
	efficiency_map_rotSCSN_81S.png \
	energyPS_all_afid.png \
	energy_PS_bins_pref_SCSN_81F1.png \
	energy_PS_bins_pref_log_SCSN_81F1.png \
	energy_PS_EJ_200.png \
	energy_PS_EJ_260.png \
	energy_PS_EJ_260_2P.png \
	energy_PS_SCSN_81F1.png \
	energy_PS_SCSN_81F2.png \
	energy_PS_SCSN_81F3.png \
	energy_PS_SCSN_81F4.png \
	energy_PS_SCSN_81S.png \
	canvas_energy_tree_EJ_200_0_1_1_lin.png \
	canvas_energy_tree_EJ_200_0_1_1_log.png \
	canvas_energy_tree_SCSN_81F1_0_1_1_lin.png \
	canvas_energy_tree_SCSN_81F1_0_1_1_log.png \
	canvas_gaussfit_lin_EJ_200.png \
	canvas_gaussfit_lin_SCSN_81F1.png \
	canvas_gaussfit_log_EJ_200.png \
	canvas_gaussfit_log_SCSN_81F1.png \
	ts.png \
	tsF.png

ADDFIGS = align_0.png \
	align_1.png \
	energyPS_all_ped_log.png \
	canvas_energy_tree_EJ_260_0_1_1_lin.png \
	canvas_energy_tree_EJ_260_0_1_1_log.png \
	canvas_energy_tree_EJ_260_2P_0_1_1_lin.png \
	canvas_energy_tree_EJ_260_2P_0_1_1_log.png \
	canvas_energy_tree_SCSN_81F1_1_1_1_lin.png \
	canvas_energy_tree_SCSN_81F1_1_1_1_log.png \
	canvas_energy_tree_SCSN_81F2_0_1_1_lin.png \
	canvas_energy_tree_SCSN_81F2_0_1_1_log.png \
	canvas_energy_tree_SCSN_81F3_0_1_1_lin.png \
	canvas_energy_tree_SCSN_81F3_0_1_1_log.png \
	canvas_energy_tree_SCSN_81F4_0_1_1_lin.png \
	canvas_energy_tree_SCSN_81F4_0_1_1_log.png \
	canvas_energy_tree_SCSN_81S_0_1_1_lin.png \
	canvas_energy_tree_SCSN_81S_0_1_1_log.png \
	canvas_gaussfit_lin_EJ_260.png \
	canvas_gaussfit_lin_EJ_260_2P.png \
	canvas_gaussfit_lin_SCSN_81F2.png \
	canvas_gaussfit_lin_SCSN_81F3.png \
	canvas_gaussfit_lin_SCSN_81F4.png \
	canvas_gaussfit_lin_SCSN_81S.png \
	canvas_gaussfit_log_EJ_260.png \
	canvas_gaussfit_log_EJ_260_2P.png \
	canvas_gaussfit_log_SCSN_81F2.png \
	canvas_gaussfit_log_SCSN_81F3.png \
	canvas_gaussfit_log_SCSN_81F4.png \
	canvas_gaussfit_log_SCSN_81S.png

$(FIGS):
	@if find . -name $@ | egrep '.*' >& /dev/null; then \
	convert `find . -name $@` -transparent white `find . -name $@`; \
	tar --transform "s|`find . -name $@ -printf '%h\n'`|figs|" \
	-uf dn-18-007_figs.tar `find . -name $@`; \
	else echo "File $@ not found!"; fi

$(ADDFIGS):
	@if find . -name $@ | egrep '.*' >& /dev/null; then \
	convert `find . -name $@` -transparent white `find . -name $@`; \
	tar --transform "s|`find . -name $@ -printf '%h\n'`|addfigs|" \
	-uf dn-18-007_figs.tar `find . -name $@`; \
	else echo "File $@ not found!"; fi

packfigs: $(FIGS) $(ADDFIGS)

cleanall: fitclean libclean plotclean
	@rm -f yield_results.txt dn-18-007_figs.tar fiduciality_test_?.png
