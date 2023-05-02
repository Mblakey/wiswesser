# WLN Parser out of Babel Build

Wisswesser WLN <---> SMILES parser

## Requirements

This project is built with CMake version > 3.00. <br>
git, cmake, make and a c++ compiler are all essential. <br>

## Build

run `build.sh` from the project directory, this will clone and build openbabel as well as linking
the library to the parser in cmake. 


## Usage

Command line utility **read-wln** should be created in the build directory of Parser. wisswesser can either take a sequence (single quote escaped) from the command line. e.g 'L6TJ'

```
read-wln <options> 'string'
```

#### flags

`-d` - enable all debugging logs to stderr
`-w` - dump the wln graph to a dot file in the build directory, this can be seen using the following commands

```
must be in build directory and have graphviz installed
dot -Tsvg wln-graph.dot -o wln-graph.svg && open wln-graph.svg
```

### Output
By default conversion will be piped to stdout, if converting files, its recommended piping the output to a file with '>'. 

