grammar WLN;

token:            ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
                  'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
atom_char:        ('B'|'C'|'K'|'M'|'N'|'O'|'P'|'S'|'V'|'W'|'X'|'Y');
terminator:       ('I'|'H'|'F'|'G'|'E'|'Q'|'Z');
digit:            ('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space:            ' ';
alkyl_chain:      digit+; 
element:          '-' token token '-'; 
hyper_valent:     '-' ('S'|'P'|'I'|'F'|'G'|'E') '-'; 
benzene:          'R'; 
locant:           token '-'*'&'*;
ring_size:        digit|'-'digit+'-'; 
fuse_bonds:       '/' locant locant; 
closures:         '&'+;
poly_subrings:    ring_size* (space locant ring_size+)*;
multi_subrings:   ring_size* fuse_bonds* (space locant ring_size+ fuse_bonds*)*;
multi_assign:     digit locant+;
multi_size:       locant '-'?; 
hetero_atoms:     (atom_char|element|'U')* (space locant (atom_char|element|'U')+)*;
multi_hetero:     locant (atom_char|element)+;
aromaticity:      ('T'|'&')*;
bridging:         (space locant)*; 
poly_cycle:       ('L'|'T') poly_subrings bridging hetero_atoms aromaticity 'J'; 
multi_cycle:      ('L'|'T') multi_subrings bridging space multi_assign space multi_size space? (multi_hetero space?)* aromaticity 'J'; 
cycle:            poly_cycle | multi_cycle; 

// each branch rule needs doubling without the branch start, splitting out the recursion here will be a challenge.
branch:           branch ('S'|element|hyper_valent) (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?  |
                  ('S'|element|hyper_valent) (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? (branch ('&'|terminator))? branch?         |
                  
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
                  ('U'|'U''U') branch                                                                                                                       | // unsaturation
                  
                  branch benzene (space locant (branch|terminator) )* closures*                                                                                           |
                  benzene (space locant (branch|terminator) )* closures*                                                                                                  |
                  (atom_char|alkyl_chain|element|hyper_valent|benzene) closures*;                                                                            
                                 

wln_branch:       terminator | (terminator? branch terminator?);  

charge:           space '&' digit+ '/' digit+;
ion:              space '&' wln_branch;

post_notation:    charge|ion; 

major_cycle:      cycle closures* (  (space locant wln_branch ('-' space locant major_cycle)? )|(space locant '-'? space locant major_cycle)|(space locant))*;  // handles spiro in the locant definition

macro_cycle:      ('L'|'T') '-' major_cycle space locant ring_size 'J'; 

read:             (wln_branch|major_cycle|macro_cycle) (post_notation)* EOF; 