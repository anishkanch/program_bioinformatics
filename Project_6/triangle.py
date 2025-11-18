#!/usr/bin/env python3


#usage should be triangle.py <character> <height>
#prints the character in a triangle format

import sys

#function 1 (pyramid) will handle calculation and printing
def pyramid(character: str, height: int) -> None:
    """
    This function will print a triangle based on the character and length provided. It takes into account even and odd numbers. 
    Even numbers have a repeted middle number peak while odds have a single peak.
    """
    middle = (height+1)//2

    #ascending, going up the pyramid
    asc = middle+1
    for count in range(1, asc):
        print(character * count)

    #descending, going down the pyramid
    if height % 2 == 0:
        dec = middle
    else:
        dec = middle - 1
    for count_2 in range(dec, 0, -1):
        print(character * count_2)


#function 2 (commandline) will pass commandline as arguments
def commandline(char: str, hei: str) -> None:
    """
    provides the argument lines in to the pyramid function
    """
    height = int(hei)
    pyramid(char, height)

#non function line
commandline(sys.argv[1], sys.argv[2])
