grammar locants; 

token: ('A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|'N'|
        'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z');
        
digit:('0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9');
space: ' ';

wln:  (token|digit|'&'|'-')+;
benzene: 'R';
benz_locant: ('A'|'B'|'C'|'D'|'E'|'F');


multisubstituted: (wln)? benzene (space benz_locant wln)*? EOF;





