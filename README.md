# WLN Parser External Build

Wiswesser Line Notation Parser, covert WLN to smiles, inchi, mol files etc
This is the development for the WLN parser that will be added directly to OpenBabel, 
on the release of this project, please see my local branch of babel for pre-release versions. 

[release](./notes/release.md)

## Requirements

**git**, **cmake**, **make** and a **c++ compiler** are all essential. <br>
**graphviz** is an optional install to view wln graphs (not needed for build). 

## Build

run `./build.sh` from the project directory, this will clone and build openbabel as well as linking
the library to the parser in cmake. 

A prompt will then ask you to run `source ./link.sh` which sets the babel build directory in shell, this tends to be needed for linux, but if openbabel is already installed via brew on macos, you can safely ignore. 

## Usage

Command line utility `readwln` should be created in the build directory of Parser. This can either take a sequence (single quote escaped) from the command line. e.g 'L6TJ'

```
readwln <options> -s 'string'
```

### flags

`-c` - enable run-time string correction, this allows correction of branching symbols and spaces where a valid WLN string can be seen <br>
`-d` - enable all debugging logs to stderr<br>
`-h` - display the help menu <br>
`-o` - choose output format for string, options are `-osmi`, `-oinchi` and `-ocan` following OpenBabels format conventions <br>
`-w` - dump the wln graph to a dot file in the build directory, this can be seen using the following commands <br>

```
dot -Tsvg wln-graph.dot -o wln-graph.svg && open wln-graph.svg
```

### Output

By default conversion will be piped to stdout, if converting files, its recommended piping the output to a file with '>'. All logging and debug information including fatal messages on unsuccesful conversion are sent to stderr. Babel Library information is also sent to stderr. 


### Unit Testing

A text file for the WLN strings contained in Elbert G. Smiths rule book are contained in data. To run the unit test, run `./smith_test.sh` contained in the src/test directory. This will give a score of successful conversions (100 is NOT expected), some incorrect strings are present. 

Data from Chembl, Chemspider and PubChem are also contained and have their corresponding unit test shell files in the `test` directory. These are run like the smith test, due to incorrect OCR it is not expected to be at 100%.

Another text file containing all the english words is included, `./english_test.sh` will parse the reader over every capatilised english word and check for seg faults, this also returns the matches to the `data` directory for word challenges. 

For Chembl, Chemspider and PubChem, run the equivilent shell files. 



