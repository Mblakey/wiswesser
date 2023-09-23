grammar WLN;

token:          ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
                'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
atom_char:      ('B'|'C'|'E'|'F'|'G'|'H'|'I'|'K'|'M'|'N'|
                'O'|'P'|'Q'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
terminators:    ('I'|'H'|'F'|'G'|'E'|'Q'|'Z');
digit:          ('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space:          ' ';
alkyl_chain:    digit+; 
element:        '-' token token '-'; 
benzene:        'R'; 
locant:         token '&'*'-'*;
ring_subset:    digit* (space locant digit+)*;
hetero_atoms:   (atom_char|element+)* (space locant (atom_char|element+))*;
aromaticity:    ('T'|'&')*;
bridging:       (space locant)*;  
poly_cycle:     ('L'|'T') ring_subset bridging hetero_atoms aromaticity 'J'; 
multi_cycle:    ('L'|'T') ring_subset bridging space digit locant+ space locant (space hetero_atoms)? aromaticity 'J'; 
cycle:          poly_cycle | multi_cycle; 

branch:         branch 'S' (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? branch?  |
                branch 'P' (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? branch?                              |
                
                branch ('X'|'K') (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&' branch)?    | 
                branch ('Y') (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&' branch)?                                           | // allows methyl contraction
                
                branch ('N'|'B') (branch ('&'|terminators))? branch?                                                                                | // no contraction allowed
                branch ('M'|'O'|'C'|'V') branch?                                                                                                    |
                branch element (branch ('&'|terminators)?)*                                                                                         |
                branch ('U'|'U''U') branch                                                                                                          | // unsaturation
                branch benzene (space locant branch)*                                                                                               | // benzene is the only allowed inline with special rules
                (atom_char|alkyl_chain|element|benzene)+;

ion:            space '&' branch;

major_cycle:    cycle (space locant branch)*; 

read:           (branch|major_cycle) (ion)* EOF; 