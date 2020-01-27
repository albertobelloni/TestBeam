SHELL := /bin/sh # sh is the default
.SILENT: clean

PHONY: clean

clean:
	mkdir -p results
	if [ -a roofit_e* ]; then mv roofit_e* results/; \
	else echo "No roofit files"; fi
	if [ -a canvas_e* ]; then mv canvas_e* results/; \
	else echo "No canvas files"; fi
	if [ -a canvas_gaussfit_* ]; then mv canvas_gaussfit* results/; \
	else echo "No canvas Gauss fit files"; fi
