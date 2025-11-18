class Read:
    """
    Stores each field of the SAM entry in named attributes 
    """
    
#Write an __init__ method to construct an instance of your Read class from a (non-header) line in a SAM file
    def __init__(self, sam_line):
        fields = sam_line.strip().split("\t")
        
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = fields[7]
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]

#Use the "flag" field of the SAM entry to return booleans corresponding to the following mapping attributes
#1 if it is true or a 0 if it is false
#https://stackoverflow.com/questions/17330160/how-does-the-property-decorator-work-in-python, to combat error of "Your Read.is_forward contains <bound method>""
        
    #1. is_mapped - 4 indicates not mapped. So false means true
    @property
    def is_mapped(self):
        return(self.flag & 4) == 0
        
    #2. is_forward - 16 indicates mapped in the reverse direction. So false means true
    @property
    def is_forward(self):
        return(self.flag & 16) == 0
    
    #3. is_reverse - 16 again, but opposite
    @property
    def is_reverse(self):
       return(self.flag & 16) != 0
    
    #4. is_primary (is this SAM entry the primary mapping of this read?)
    #256 - not the best mapping should be false (AKA, best mapping is true)
    #2048 - larger part mapped elsewhere should be false
    @property
    def is_primary(self):
        return(self.flag & 256) == 0 and (self.flag & 2048) == 0

#Write a base_at_pos method to take a 1-base position in the reference sequence and return the base mapped to specific position in the reference.
    def base_at_pos(self, pos: int) -> str:
        #for empty cigar or position out of range
        if self.cigar == "*" or pos < self.pos:
            return ""
        
        reference = self.pos
        reading = 0
        num = ''

        #obtain number and letter, for example if it is 3M it will isolate 3 and M
        #unsure if usage of parsar is allowed, studied parsar and create a parsar code instead of using library
        #https://stackoverflow.com/questions/607760/python-parsing 
        for unit in self.cigar:
            if unit.isdigit():
                num = num + unit
            else:
                number = int(num)
                letter = unit
                num = ''
        
                #return perfect matches if they fall inside matching region
                if letter == "M":
                    if reference + number > pos:
                        return self.seq[reading + (pos - reference)]
                    reading = reading + number
                    reference = reference + number
                
                #add  on until matches are found
                elif letter == "I":
                    if reference - 1 == pos or reference == pos:
                        return self.seq[reading: (reading + number)]
                    reading = reading + number
                elif letter == "D":
                    if reference + number > pos:
                        return ""
                    reference = reference + number
                elif letter == "S":
                    reading = reading + number 
                elif letter == "H":
                    continue
        
        #if no matches, return empty string
        return ""

#Write a qual_at_pos method to return the quality score of a base mapped to a specific position in the reference
#will be the same as base_at_pos, but replace self.seq with self.qual
    def qual_at_pos(self, pos: int) -> str:
        #for empty cigar or position out of range
        if self.cigar == "*" or pos < self.pos:
            return ""
        
        reference = self.pos
        reading = 0
        num = ''

        #obtain number and letter, for example if it is 3M it will isolate 3 and M
        for unit in self.cigar:
            if unit.isdigit():
                num = num + unit
            else:
                number = int(num)
                letter = unit
                num = ''
        
                #return quality if they are perfect matches and if they fall inside matching region
                if letter == "M":
                    if reference + number > pos:
                        return self.qual[reading + (pos - reference)]
                    reading = reading + number
                    reference = reference + number
                
                #add  on until matches are found
                elif letter == "I":
                    reading = reading + number
                elif letter == "D":
                    reference = reference + number
                elif letter == "S":
                    reading = reading + number 
                elif letter == "H":
                    continue
        
        #if no matches, return empty string
        return ""

#Write a mapped_seq method which returns the mapped portion of the read that corresponds to the reference sequenced
    def mapped_seq(self) -> str:
        #for empty cigar
        if self.cigar == "*":
            return ""
        
        reference = self.pos
        reading = 0
        num = ''
        mapped_portion = ''

        #obtain number and letter, for example if it is 3M it will isolate 3 and M
        for unit in self.cigar:
            if unit.isdigit():
                num = num + unit
            else:
                number = int(num)
                letter = unit
                num = ''
        
                #build a list based on matches, insertions, deletions... ignore S and H
                if letter == "M":
                    mapped_portion = mapped_portion + self.seq[reading:(reading+number)]
                    reading = reading + number
                    reference = reference + number
                elif letter == "I":
                    mapped_portion = mapped_portion + self.seq[reading:(reading+number)]
                    reading = reading + number
                elif letter == "D": #per directions, include dash if deletion
                    mapped_portion = mapped_portion + ('-' * number)
                    reference = reference + number
                elif letter == "S":
                    reading = reading + number 
                elif letter == "H":
                    continue
        
        #return final mapped poriton
        return mapped_portion
