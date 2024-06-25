# `findTDE` upcoming work
* remove dependence on slurm submission script and make more transferrable
-> how to do: need most of script to operate in typical directory and exclusively perform VASP calculations outside
    * script itself may be fine, just replacing slurm submission with vasp execution
    * maybe use flags in pbs script to determine when to perform elsewhere: execute find_tde until specific comment line, perform a command, then continue original script
    * could set up multiple functions then on script call, execute certain functions
    * option to specify separate directory for executing vasp?
* add in VASP ML utility

## multi_tde.py
* move examples to an example doc
* maybe change from screen? dunno

## tde_analysis.py
* remove built in requirement for ga 34/etc., require input for specific things
* 1 atom displacement at a time
* maintain plotting defaults, optional user input
* examples with lines for GaN displacements to copy-paste for auto analysis

## defect_analysis.py
* remove built in dphi dtheta for averaging TDEs
* double check defect type checking
* rework defect location analysis