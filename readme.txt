This is the top level readme for the NEURON models associated with the
paper:

Hsu CL, Zhao X, Milstein AD, Spruston N (2018) Persistent sodium
current mediates the steep voltage dependence of spatial coding in
hippocampal pyramidal neurons. Neuron 99, 1–16

These models were based on previously published models (as described
in the paper). The simulations serve as independent tests for the
proposed mechanism underlying the synaptic amplification observed in
our experiments.

These models require that the NEURON simulation environment is
installed, which is freely available from http://neuron.yale.edu

To run the simple model see the instructions in the readme in the
SimpleModel folder, and for the full model (which is also available at
the github repo https://github.com /neurosutras/nap) see the readme in
the FullModel folder.

These models were contributed by Ching-Lung Hsu
(hsuc@janelia.hhmi.org) and Aaron Milstein (aaronmil@stanford.edu).

20180615 An update from Aaron Milstein for the FullModel folder where
the three python scripts for generating figures take in the required
.hdf5 file paths as command line arguments, e.g.:

python plot_nap_EPSC_amplification data/20180220_nap_EPSC_amplification_soma.hdf5  # generates plots similar to Figure S5
python plot_nap_EPSP_amplification data/07162017_nap_amplification_DC_soma_stim_trunk.hdf5 data/07162017_nap_amplification_DC_dend_stim_trunk.hdf5  # generates plots similar to Figure 5C and S4
python plot_nap_EPSP_amplification_IO data/08112017_nap_amplification_IO_DC_soma_stim_trunk.hdf5  # generates plots similar to Figure 6B
