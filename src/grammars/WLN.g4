grammar WLN;

token:            ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
                  'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
atom_char:        ('B'|'C'|'E'|'F'|'G'|'H'|'I'|'K'|'M'|'N'|
                  'O'|'P'|'Q'|'S'|'T'|'U'|'V'|'W'|'X'|'Y');
terminator:      ('I'|'H'|'F'|'G'|'E'|'Q'|'Z');
digit:            ('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space:            ' ';
alkyl_chain:      digit+; 
element:          '-' token token '-'; 
benzene:          'R'; 
locant:           token '-'*'&'*;
ring_size:        digit|'-'digit+'-'; 
fuse_bonds:       '/' locant locant; 
closures:         '&'+;
poly_subrings:    ring_size* (space locant ring_size+)*;
multi_subrings:   ring_size* fuse_bonds* (space locant ring_size+ fuse_bonds*)*;
multi_assign:     digit locant+;
multi_size:       locant; 
hetero_atoms:     (atom_char|element)* (space locant (atom_char|element)+)*;
multi_hetero:     locant (atom_char|element)+;
aromaticity:      ('T'|'&')*;
bridging:         (space locant)*; 
poly_cycle:       ('L'|'T') poly_subrings bridging hetero_atoms aromaticity 'J'; 
multi_cycle:      ('L'|'T') multi_subrings bridging space multi_assign space multi_size space (multi_hetero space?)* '-'? aromaticity 'J'; 
cycle:            poly_cycle | multi_cycle; 

// each branch rule needs doubling without the branch start, splitting out the recursion here will be a challenge.
branch:           branch ('S'|element) (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?  |
                  ('S'|element) (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?         |
                  
                  branch 'P' (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?                                       |
                  'P' (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?                                              |
                
                  branch ('X'|'K') (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&' branch)?             | 
                  ('X'|'K') (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&' branch)?                    | 
                  
                  branch ('Y') (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&' branch)?                                                   | // allows methyl contraction
                  ('Y') (branch ('&'|terminator)? | '&')? (branch ('&'|terminator)? | '&' branch)?                                                          | // allows methyl contraction
                
                  branch ('N'|'B') (branch ('&'|terminator))? (branch ('&'|terminator))? branch?                                                            | // no contraction allowed
                  ('N'|'B') (branch ('&'|terminator))? (branch ('&'|terminator))? branch?                                                                   | // no contraction allowed
                  
                  branch ('M'|'O'|'C'|'V') branch?                                                                                                          |                                                                                        
                  ('M'|'O'|'C'|'V') branch?                                                                                                                 |                                                                                        

                  branch ('U'|'U''U') branch                                                                                                                | // unsaturation
                  
                  branch benzene (space locant branch)* closures*                                                                                           |
                  benzene (space locant branch)* closures*                                                                                                  |
                  (atom_char|alkyl_chain|element|benzene) closures*;                                                                                           

wln_branch:       terminator? branch terminator?;  
ion:              space '&' wln_branch;

major_cycle:      cycle closures* (  (space locant wln_branch ('-' space locant major_cycle)? )|(space locant '-'? space locant major_cycle)|(space locant))*;  // handles spiro in the locant definition

read:             (wln_branch|major_cycle) (ion)* EOF; 