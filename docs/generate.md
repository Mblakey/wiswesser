# WLN Generation 

This file contains usage notes for generating compounds in WLN that trend toward a target chemical descriptor. This software uses the Finite State Machine (FSM) as a markov decision model, with Q-learning applied to a non-deep learning method to allow rapid generation and machine-learning where the black box doesn't hold us back. <br>

## Executables

Command line utility `wlngen`. Most of the arguments here are parameters for the Q-learning process, and follow a `<flag>=<value>` format. Any input files are used to seed the probabilities in the markov chain, and are highly recommended to speed up useful convergence.<br> 

The file format for the train files should be WLN strings split with the newline character, examples of which are contained in `data`. <br>

From the build directory:<br>

```
./wlngen <options> <filename(s)>
```

#### Flags 

Standard Flags: <br>
`-l=<int>` or `--length=<int>` - set length for generation (default 5) <br>

Hyperparameter Flags: <br>
`-r=<int>` or `--runs=<int>` - set learning episodes (default 5)<br>
`-e=<double>` or `--epsilon=<double>` - set epsilon hyperparameter (default 0.5) <br>
`-d=<double>` or `--decay=<double>` - set decay rate hyperparameter (default 0.005) <br>
`-a=<double>` or `--alpha=<double>` - set learning rate hyperparameter(default 0.5) <br>
 
General Flags: <br>
`-p` or `--print` - show all set hyperparameters and exit <br>
`-h` or `--help` -show this help menu and exit <br>
