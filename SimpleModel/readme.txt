This is the readme for the simple model described in the paper:

Hsu CL, Zhao X, Milstein AD, Spruston N (2018) Persistent sodium
current mediates the steep voltage dependence of spatial coding in
hippocampal pyramidal neurons. Neuron 99, 1–16

This model requires NEURON which is freely available from
http://neuron.yale.edu

Usage:

Compile the mod files with a command like "nrnivmodl" (linux/unix) or
use mknrndll (mswin or mac os x).

Run by executing any of the startX.hoc files where X=CC, InOutCurve
or VC. An example command on the linux/unix platform:

nrngui StartVC.hoc

After the simulation has completed, plots similar to Figure 5A (X=CC),
5B (X=VC) or 6A (X=InOutCurve) will be generated.
