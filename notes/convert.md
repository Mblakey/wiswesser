

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