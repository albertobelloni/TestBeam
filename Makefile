SHELL := /bin/sh # sh is the default
.SILENT: clean

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
