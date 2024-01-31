# WLN Extraction 

This file contains usage notes for extracting WLN from any parsable ascii encoded document, it uses an Finite State Machine (FSM) with a pushdown stack in order to very quickly parse documents. If you are familiar with the `grep` tool, flags follow the convection set. 


## Executables

Command line utility `wlngrep`. This either takes a filename, or an escaped single sequence if using the `-s` flag. With some exceptions, flags are kept inline with standard grep usage.

```
wlngrep <options> <filename>
```

#### Flags 

`-c` - return number of matches instead of string <br>
`-o` - print only the matched parts of line <br>
`-m` - do not minimise DFA (debugging only) <br>
`-s` - interpret `<filename>` as a string to match <br>
`-x` - return string if whole line matches <br>