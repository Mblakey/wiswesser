# The Wiswesser Line Notation (WLN) Project 

* WLN Parser               - read and write WLN to/from smiles, inchi, mol files and other chemical line notations.  <br>
* WLN FSM                  - extract chemical terms from documents, this machine uses greedy matching to return matched WLN sequences from documents.  <br>
* WLN Compresser      - compress WLN strings using markov decision processes. <br>

This is Linux and MacOS software only. <br>

Note: This project is solely created by Michael as part of his PhD work, if you are interested using the project, or find any bugs or issues, reporting them would be extremely helpful. 

<img src="./docs/intro_wln.png" width="800">

## Requirements

**git**, **cmake**, **make** and a **c++ compiler** are all essential. <br>
**graphviz** is an optional install to view wln graphs (not needed for build). 

**OpenBabel** [see repo](https://github.com/openbabel/openbabel), will be installed as an external dependency.  

## Build

Run `./bootstrap.sh` from the project directory, this will clone and build openbabel as well as linking
the library to the parser in cmake. Babel files will be installed to `external`. Building the projects places all executables into `build/`. <br>

## Project Structure 

This repository contains a broad range of functionality using WLN notation for various operations. As such, please read the individual `README.txt` files for the required area. <br>

* [WLN Conversion - Read and Write](./docs/convert.md) <br>
* [WLN Compression - Lossless Compress/Decompress](./docs/compress.md)<br>
* [WLN Extract - FSM wlngrep](./docs/extract.md) <br>


## Unit Testing

All unit tests are contained in the `/test` directory. <br>
These include: 
1. `compare.sh`
2. `reading.sh`
3. `writing.sh`
4. `file.sh`

Unit tests 1-3 operate on the data files in `\data`. For comparsions agaisnt the old parser in OpenBabel select 1, for reading count tests run 2, writing round trip tests 3. To parse a file of WLN strings, `file.sh` will attempt conversions on every line.


