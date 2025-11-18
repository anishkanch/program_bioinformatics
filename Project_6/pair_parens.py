#!/usr/bin/env python3

#usage should be pair_parens.py <test string>

import sys

#function to return paired on unpaired
def parenthesis (input: str) -> None:
    """
    Prints PAIRED if all open parenthesis are followed by closed parentheses or paired with one later.
    Prints NOT PAIRED if there are lone parentheses, open or closed.
    """

    count = 0
    for i in input:
        if i == "(":
            count = count+1
        elif i == ")":
            count = count-1
            if count < 0:
                print("NOT PAIRED")
                return
    
    if count == 0:
        print("PAIRED")
    else:
        print("NOT PAIRED")

#function to send commandline into parenthesis
def commandline(argument: str) -> None:
    """
    This function will take the parenthesis string and provide it to the parenthesis function previously listed.
    """
    parenthesis(argument)

#no function, read in argument
commandline(sys.argv[1])
