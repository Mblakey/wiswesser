grammar WLN;

token:  ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
        'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
atom_char:  (  'B'|'C'|'E'|'F'|'G'|'H'|'I'|'K'|'M'|'N'|
              'O'|'P'|'Q'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
terminators: ('I'|'H'|'F'|'G'|'E'|'Q'|'Z');

digit:('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space: ' ';
dash: '-';
slash:'/';

carbon_chain: digit+; 
element: dash token token dash; 

branch:   branch 'S' (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? |
          branch 'P' (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? (branch ('&'|terminators))? |
          branch ('X'|'K') (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&' branch)? | 
          branch ('Y') (branch ('&'|terminators)? | '&')? (branch ('&'|terminators)? | '&' branch)? | // allows methyl contraction
          branch ('N'|'B') (branch ('&'|terminators))? (branch ('&'|terminators))? | // no contraction allowed
          branch ('M'|'O'|'C'|'V') branch?  |
          branch element (branch '&'? | '&')* |
          branch 'U'+ branch | // unsaturation
          (atom_char|carbon_chain|element)+ ;


ion: space '&' branch;


ring: EOF;

read: (branch|ring) (ion)* EOF; 