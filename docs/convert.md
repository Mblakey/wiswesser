# Converting between WLN and CLN Formats

This file contains usage notes for coverting between WLN and more modern chemical line notations, also, experimentally, a modern version of WLN used for the compression and generation projects. 

## Executables


### `readwln`

`readwln` - This takes a WLN sequence (single quote escaped) from the command line. e.g 'L6TJ', and returns the desired output format.<br> 

From the build directory:<br>

```
./readwln <options> -o<format> 'string'
```

#### Flags

`-h` - display the help menu <br>
`-o` - choose output format for string, options are `-osmi`, `-oinchi`, `-okey` (inchikey)and `-ocan` following OpenBabels format conventions <br>
`--old` - use nextmoves old wln parser (lower coverage, much faster)<br>


### `writewln`

`writewln` - This takes an input sequence (single quote escaped) from the command line. e.g 'c1ccccc1', and returns the corresponding WLN string.<br> 

From the build directory:<br>

```
./writewln <options> -i<format> 'string'
```

#### Flags 

`-h` - display the help menu <br>
`-i` - choose input format for string, options are `-ismi`, `-iinchi` and `-ican` following OpenBabels format conventions <br>
`-m` - generate modern WLN notation (experimental)


## Wiswesser Conversion Release Notes

The following are sections from Elbert G. Smiths rule book that were used to create the wln reader. Note that not all chapters are listed here, only the ones where compound types were introduced.

Please note that the "MANTRAP" rules, are not officialy given in either volume of the offical Wiswesser manuals, as such, implementation is tricky at best. For this parser, they will not be supported. 

| Rule | Read | Write |
| ---- | ---- | ---- |
|Unbranched and Branched Chains | :heavy_check_mark: | :heavy_check_mark: | 
|Systematic Contractions | :heavy_check_mark: | :heavy_check_mark: |
|Organic Salts | :heavy_check_mark: | :heavy_check_mark: |
|Benzene Derivatives | :heavy_check_mark: | :heavy_check_mark: |
|Multisubstituted Benzene Rings | :heavy_check_mark: | :heavy_check_mark: |
|Benzene Rings in Branching Chains | :heavy_check_mark: | :heavy_check_mark: |
|Monocyclic Rings | :heavy_check_mark: | :heavy_check_mark: |
|Bicyclic Rings | :heavy_check_mark: | :heavy_check_mark: |
|Polycyclic Rings | :heavy_check_mark: | :heavy_check_mark: |
|Perifused Rings | :heavy_check_mark: | :heavy_check_mark: |
|Chains of Rings other than Benzene | :heavy_check_mark: | :heavy_check_mark: |
|Sprio Rings | :heavy_check_mark: | :heavy_check_mark: |
|Bicyclic Bridged Rings |:heavy_check_mark: | :heavy_check_mark: |
|Rings with Pseudo Bridges | :heavy_check_mark: | :heavy_check_mark: |  
|Ring Structures with Crossed Bonds and Unbranched Bridges | :heavy_check_mark: | :heavy_check_mark: |
|Rings of Rings Contraction | :heavy_check_mark: |  | 
|Metallocenes and Catanenes | :heavy_check_mark: | :heavy_check_mark: |
|Chelete Compounds | :heavy_check_mark: | :heavy_check_mark: |
|Ionic Charges, Free Radicals and Isotopes | :heavy_check_mark: | :heavy_check_mark: |
|Multipliers | | |
|Ring Contractions and Multipliers | | |
|All Special Problems Rules | | |
