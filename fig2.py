### Analysis of DG network data ###
# This Python code creates a scatter plot of output vs input sim scores.
# Enter the idname

# ModelDB file along with publication:
# Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
# http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

# modified and augmented by
# Man Yi Yim / 2015
# Alexander Hanuschkin / 2011

idname = "-pp10-gaba1-kir1-st30"
execfile('plot_DG_all.py')
execfile('GCinput.py')
execfile('inout_pattern.py')