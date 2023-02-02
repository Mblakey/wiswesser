# WLN Parser out of Babel Build

Wisswesser WLN <---> SMILES parser

## Requirements

This project is built with CMake version > 3.00. <br>
OpenBabel must be built and linked in the CMakeLists.txt. 

## Build

**Openbabel** must be dynamically linked in the CMakeLists.txt, static build for obabel is incompatible with maeparser at current time. CMakeLists.txt expects openbabel to be in the CodeBase folder and compiled in build. 

```
// -- babel set up --

cd CodeBase
git clone https://github.com/openbabel/openbabel.git
cd openbabel 
mkdir build
cd build
cmake ..
make -j10

// -- parser --

cd CodeBase
cd Parser
mkdir build
cd build
cmake ..
make
```


## Usage

Command line utility **wisswesser** should be created in the build directory of Parser. wisswesser can either take a input file or a sequence single quote escaped from the command line. e.g 'L6TJ' or 'c1ccccc1'.

### Modes
Mandatory arguments include \<input\> followed by **--i\<format\>** to specify input format and **--o\<format\>** to specify output conversion.  

| mode | string type|
| --- | ------ |
| smi | smiles | 
| wln | wln    | 
| inc| inchi  |
| inckey| inchi key |


###Output
By default conversion will be piped to stdout, if converting files, its recommended piping the output to a file with '>'. 

###Example 

```
wisswesser input_wln.txt --iwln --osmi > converted_smiles.txt
```