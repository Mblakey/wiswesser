grammar wln_fsm;

// Removes all semantics, only syntactics 

token:            ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
                  'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');

contract_char:    ('X'|'K'|'Y'); // methly contractions allowed

branch_char:      ('B'|'N'|'P'|'S'); // allowed branches

linear_char:      ('C'|'O'|'V'|'M'); // not allowed branches, restrict dioxo sepereately

standard_set:     (contract_char|branch_char|linear_char); 

terminator:       ('I'|'H'|'F'|'G'|'E'|'Q'|'Z'); // terminates branch
digit:            ('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space:            ' ';
alkyl_chain:      digit+; 
element:          '-' token token '-'; 
hyper_valent:     '-' ('S'|'P'|'I'|'F'|'G'|'E') '-'; 
benzene:          'R'; 

allowed_dioxo_start:  'W' (contract_char|branch_char|element|hyper_valent); 

/* OG 
// restrict 'W' here
branch:   (contract_char|branch_char|element|hyper_valent) 'W'? branch? |
          (linear_char|alkyl_chain) branch? |
          terminator branch?  |
          (standard_set|alkyl_chain|terminator) '&'+ branch? | // will have two eliminate the recursion fully
          ('U'|'U''U') branch |
          benzene (benzene+|'&'|EOF)
          ;
*/     

branch:   (contract_char|branch_char|element|hyper_valent) 'W'? branch? |
          (linear_char|alkyl_chain) branch? |
          (standard_set|alkyl_chain|terminator) '&'+ branch? | // will have two eliminate the recursion fully
          terminator branch?  |
          ('U'|'U''U') branch |
          benzene (benzene+|'&'|EOF)
          ;

substitued_benzene: benzene (space locant branch)*;

wln_branch:       allowed_dioxo_start? (branch|substitued_benzene)*;  



/* ############## CYCLES ############## */

locant:           token '-'*'&'*;
ring_size:        digit|'-'digit+'-'; 
fuse_bonds:       '/' locant locant; 
closures:         '&'+;
poly_subrings:    ring_size* (space locant ring_size+)*;
multi_subrings:   ring_size* fuse_bonds* (space locant ring_size+ fuse_bonds*)*;
multi_assign:     digit locant+;
multi_size:       locant '-'?; 
hetero_atoms:     (standard_set|element|'U')* (space locant (standard_set|element|'U')+)*;
multi_hetero:     locant (standard_set|element)+;
aromaticity:      ('T'|'&')*;
bridging:         (space locant)*; 
poly_cycle:       ('L'|'T') poly_subrings bridging hetero_atoms aromaticity 'J'; 
multi_cycle:      ('L'|'T') multi_subrings bridging space multi_assign space multi_size space? (multi_hetero space?)* aromaticity 'J'; 
cycle:            poly_cycle | multi_cycle; 


/* ############## BUILD ############## */

charge:           space '&' digit+ '/' digit+;
ion:              space '&' wln_branch;

post_notation:    charge|ion; 


major_cycle:      cycle closures* (  (space locant)|(space locant wln_branch ('-' space locant major_cycle)? )|(space locant '-'? space locant major_cycle))*;  // handles spiro in the locant definition


macro_cycle:      ('L'|'T') '-' major_cycle space locant ring_size 'J'; 

read:             (wln_branch|major_cycle|macro_cycle) (post_notation)* EOF; 