# WLN Parser out of Babel Build

Wisswesser Line Notation Parser, covert WLN to smiles, inchi, mol files etc
This is the developement for the WLN parser that will be added directly to OpenBabel, 
on the release of this project, please see my local branch of babel for pre-release versions. 

## Requirements

**git**, **cmake**, **make** and a **c++ compiler** are all essential. <br>
**graphviz** is an optional install to view wln graphs. 

## Build

run `build.sh` from the project directory, this will clone and build openbabel as well as linking
the library to the parser in cmake. 


## Usage

Command line utility `read-wln` should be created in the build directory of Parser. wisswesser can either take a sequence (single quote escaped) from the command line. e.g 'L6TJ'

```
read-wln <options> 'string'
```

### flags

`-d` - enable all debugging logs to stderr<br>
`-w` - dump the wln graph to a dot file in the build directory, this can be seen using the following commands

```
dot -Tsvg wln-graph.dot -o wln-graph.svg && open wln-graph.svg
```

### Output
By default conversion will be piped to stdout, if converting files, its recommended piping the output to a file with '>'. 

