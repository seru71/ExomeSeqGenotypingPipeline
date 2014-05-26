
"""
Miscellaneous utility functions

author: Pawel Sztromwasser (pawel.sztromwasser at k2.uib.no)
"""

def quote_aware_split(string, delim=',', quote='"'):
    """ Split outside of quotes (i.e. ignore delimiters within quotes."""
    out = []
    s = ''
    open_quote=False
    for c in string:
        if c == quote: 
            open_quote = not open_quote
        if c == delim and not open_quote:
            out += [s]
            s = ''
        else: 
            s += c
    return out + [s]


def parenthesis_aware_split(string, delim=',', open_par='(', close_par=')'):
    """ Split outside of parenthesis (i.e. ignore delimiters within parenthesis."""
    out = []
    s = ''
    open_parenthesis=0
    for c in string:
        if c == open_par: 
            open_parenthesis+=1
        if c == close_par and open_parenthesis > 0:
            open_parenthesis-=1
        if c == delim and open_parenthesis==0:
            out += [s]
            s = ''
        else: 
            s += c
    return out + [s]

