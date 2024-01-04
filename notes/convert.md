

# Converting between WLN and CLN Formats

This file contains usage notes for coverting between WLN and more modern chemical line notations, also, experimentally, a modern version of WLN used for the compression and generation projects. 

## Executables

`readwln` - This takes a WLN sequence (single quote escaped) from the command line. e.g 'L6TJ', and returns the desired output format. 

```
readwln <options> -o<format> -s 'string'
```

#### Flags

`-c` - enable run-time string correction, this allows correction of branching symbols and spaces where a valid WLN string can be seen <br>
`-d` - enable all debugging logs to stderr<br>
`-h` - display the help menu <br>
`-o` - choose output format for string, options are `-osmi`, `-oinchi` and `-ocan` following OpenBabels format conventions <br>
<br>

`writewln` - This takes an input sequence (single quote escaped) from the command line. e.g 'c1ccccc1', and returns the corresponding WLN string. 


```
writewln <options> -i<format> -s 'string'
```

#### Flags 

`-d` - enable all debugging logs to stderr<br>
`-h` - display the help menu <br>
`-i` - choose input format for string, options are `-ismi`, `-iinchi` and `-ican` following OpenBabels format conventions <br>
`-m` - generate modern WLN notation (experimental)


## Wiswesser Conversion Release Notes

The following are sections from Elbert G. Smiths rule book that were used to create the wln reader. Note that not all chapters are listed here, only the ones where compound types were introduced.

| Rule | Read | Write |
| ---- | ---- | ---- |
|Unbranched and Branched Chains | :heavy_check_mark: | | 
|Systematic Contractions | | |
|Organic Salts | | |
|Benzene Derivatives | | |
|Multisubstituted Benzene Rings | | |
|Benzene Rings in Branching Chains | | |
|Monocyclic Rings | | |
|Bicyclic Rings | | |
|Polycyclic Rings | | |
|Perifused Rings | | |
|Chains of Rings other than Benzene | | |
|Sprio Rings | | |
|Bicyclic Bridged Rings | | |
|Rings with Pseudo Bridges | | |  
|Ring Structures with Crossed Bonds and Unbranched Bridges | | |
|Rings of Rings Contraction | | | 
|Metallocenes and Catanenes | | |
|Chelete Compounds | | |
|Ionic Charges, Free Radicals and Isotopes | | |
|Multipliers | | |
|Ring Contractions and Multipliers | | |
|All Special Problems Rules | | |



## Rules fully supported
* Unbranched and Branched Chains
* Systematic Contractions
* Organic Salts
* Benzene Derivatives
* Multisubstituted Benzene Rings
* Benzene Rings in Branching Chains
* Monocyclic Rings
* Bicyclic Rings
* Polycyclic Rings
* Perifused Rings
* Chains of Rings other than Benzene
* Sprio Rings
* Bicyclic Bridged Rings
* Rings with Pseudo Bridges
* Ring Structures with Crossed Bonds and Unbranched Bridges
* Chelete Compounds
* Metallocenes and Catanenes 
* Ionic Charges, Free Radicals and Isotopes

