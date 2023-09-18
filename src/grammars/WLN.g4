grammar WLN;

token:  ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
        'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');

atom_set:  (  'B'|'C'|'E'|'F'|'G'|'H'|'I'|'K'|'M'|'N'|
              'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');

terminators: ('I'|'H'|'F'|'G'|'E'|'Q'|'Z');

digit:('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space: ' ';
dash: '-';
slash:'/';
element: dash token token dash; 

// for a push down automaton, track '&' so no need for strict grammar rule
// could opt repeat, but lose all branch strictness if inf allowed.

// branch recursions need to come in pairs -> starting + non starting

// split out methyl contraction rules for X,Y,K. 
branch:   ((atom_set|digit|element)('&')?)+ | // branch can be optionally terminated 
          branch ('U')* branch | // opt repeat on unsaturations
          ('X'|'K') (branch|'&')? (branch|'&')? (branch|'&')? (branch|'&')? |
          branch ('X'|'K') (branch|'&')? (branch|'&')? (branch|'&')? (branch|'&')? |
          
          ('Y'|'N'|'B') (branch|'&')? (branch|'&')? (branch|'&')? |
          branch ('Y'|'N'|'B') (branch|'&')? (branch|'&')? (branch|'&')?

          element (branch|'&')* | // let hypervalence inf on branch
          branch element (branch|'&')*; 
         

ring: EOF;

read: (branch|ring) EOF; 