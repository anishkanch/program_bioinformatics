

import re

class Read:
	def __init__(self, sam_line: str):
		(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags) = sam_line.strip().split("\t")

		# store basic properties of the read
		self.qname: str = qname
		self.flag: int = int(flag)
		self.rname: str = rname
		self.pos: str = int(pos)
		self.mapq: str = int(mapq)
		self.cigar: str = cigar
		self.rnext: str = rnext
		self.pnext: str = int(pnext)
		self.tlen: str = int(tlen)
		self.seq: str = seq
		self.qual: str = qual
		self.tags: list[str] = tags

		# score mapping properties based on flag
		self.is_mapped: bool = not bool(self.flag & 4) # 4 bit not in flag
		self.is_forward: bool = not bool(self.flag & 16)
		self.is_reverse: bool = bool(self.flag & 16)
		self.is_primary: bool = not (bool(self.flag & 256) or bool(self.flag & 2048)) # not secondary or supplemental

		# add data for mapped reads only
		self.cigar_bits: tuple[tuple[int, str]] = None
		self.mapped_len: int = None
		
		if self.is_mapped:
			self.cigar_bits = tuple([(int(n), cig) for n, cig in re.findall(r"(\d+)([A-Z])", self.cigar)])
			self.mapped_len = sum([n for n, cig in self.cigar_bits if cig in {"M", "D"}])
		

	def read_idx_at_pos(self, pos: int) -> list[None|int]:
		if not self.is_mapped:
			return []
		
		# adjust by read start
		pos -= self.pos
		if pos < 0: # If read mapped to the right of requested location
			return []

		# Check if the requested position is right of our read
		if pos >= self.mapped_len:
			return []
		
		cigar_bits = re.findall(r"(\d+)([A-Z])", self.cigar)
		# pad position based on cigar
		pad = 0
		mapped_count = 0
		for n, (size, cig_type) in enumerate(cigar_bits):
			size = int(size)
			if cig_type in {"S", "H", "I"}:
				pad += size
				continue
			if mapped_count + size >= pos+1:
				if cig_type == "M":
					# pad remaining
					pad += (pos-mapped_count)
					break
				if cig_type == "D":
					return []
			else:
				if cig_type == "M":
					mapped_count += size
					pad += size
				if cig_type == "D":
					mapped_count += size

		# Check if next bases are insertion
		if mapped_count + size == pos+1:
			if n+1 != len(cigar_bits):
				size, cig_type = cigar_bits[n+1]
				size = int(size)
				if cig_type == "I":
					return [i for i in range(pad, pad+size+1)]

		return [pad]


	def mapped_seq(self) -> str:
		if not self.is_mapped:
			return ""

		idx = 0 # track where we are in read seq
		bases = [] # list to build up over time without costly string concatenation
		for n, cig in self.cigar_bits:
			if cig == "S":
				idx += n
			elif cig == "D":
				bases += ["-"]*n
			elif cig in {"M", "I"}:
				bases += [self.seq[i] for i in range(idx, idx+n)]
				idx += n

		return "".join(bases)

	def base_at_pos(self, pos: int) -> str:
		idx = self.read_idx_at_pos(pos)
		return "".join([self.seq[i] for i in idx])

	
	def qual_at_pos(self, pos: int) -> str:
		idx = self.read_idx_at_pos(pos)
		return "".join([self.qual[i] for i in idx])
	

class SAM:
	def __init__(self):
		self.references = {} #dictionary to store reference names and lengths
		self.reads= [] #read objects
	
	@classmethod
	def from_sam(cls, sam_path):
		"""
		Load sam file, return same object
		Creates instance of Read() class
		"""
		
		sam_object = cls()

		with open(sam_path, "r") as file:
			for line in file:
				line = line.strip()

				#reference information
				if line.startswith("@SQ"):
					lin = line.split('\t')
					
					name = 'x'
					length = 0

					for li in lin:
						if li.startswith("SN:"): #names in file follow format "SN: name"
							name = li[3:]
						elif li.startswith("LN:"): #length in file follow format "LN: length"
							length = int(li[3:])
					
					if name != 'x' and length != 0:
						sam_object.references[name] = length
					
					continue
		
				#skip headers, starts with @
				if line.startswith("@"):
					continue

				#create Read() instance
				sam_object.reads.append(Read(line))
		
		return(sam_object)
	
	def reads_at_pos(self, position):
		"""
		Return a list of reads that map to a position in the reference
		We can use read_idx_at_pos from the Read class given in last weeks answer
		"""

		instances = []
		
		for reads in self.reads:
			
			#handle secondary and supplemental 
			if reads.flag & 0x100 != 0: 
				continue
			if reads.flag & 0x800 != 0:
				continue
			if not reads.is_mapped: #not mapped reads
				continue

			#position must be in reference span
			if not(reads.pos <= position < reads.pos + reads.mapped_len):
				continue
			
			index = reads.read_idx_at_pos(position) 
			
			#per previous Class, [] means deletion
			if index == []:
				instances.append(reads)
				continue

			if index is not None:
				instances.append(reads)
		
		return(instances)

	def pileup_at_pos(self, seq_name, position):
		"""
		Return a tuple containing base, qual, read from overlapping at postion, method above
		All methods we have written last week
		"""
		bases = []
		qualities = []

		overlap = self.reads_at_pos(position) #previous method

		for read in overlap:
			if read.rname.split()[0] == seq_name.split()[0]:
				base = read.base_at_pos(position) #last assignment
				quality = read.qual_at_pos(position) #last assignment

				bases.append(base)
				qualities.append(quality)
			else:
				continue
		
		return(bases, qualities)
	
	def consensus_at_pos(self, seq_name, position):
		"""
		Give a sequence name, return the most common base at a position (>50%)
		"""
		bases, qualities = self.pileup_at_pos(seq_name, position) #pass into previous function

		#unmapped positions should return empty string
		if len(bases) == 0:
			return ""
		
		#deletion vs not
		d_count = 0
		bases_list = []

		for b in bases:
			if b== "" or b== "-":
				d_count += 1
			else:
				bases_list.append(b)

		#if any base's counts are greater than 50% of the total counts, that is our consensus
		total = len(bases)
		fifty = total / 2

		if d_count > fifty:
			return ""
		if len(bases_list) == 0:
			return "N"

		#create a dictionary to story base and its corresponding count
		base_counts = {}
		#sum up the counts per base
		for base in bases_list:
			if base in base_counts:
				base_counts[base] += 1
			else:
				base_counts[base] = 1

		#find bases that are over 50%
		for base in base_counts:
			if base_counts[base] > fifty:
				return(base)
		
		return("N") #return nothing if no base is over 50%
	
	def consensus(self,seq_name):
		"""
		return majority base call from all posotions in reference (innit)
		1 - return '' if sequence name does not exist in reference
		2 - reurn '' if no reads map to it
		"""

		#1 blank return clause
		if seq_name not in self.references: #self_references is from the innit method created at the front
			return ""
		
		consensus = ''
		mapped = False

		for position in range(1, (self.references[seq_name] + 1)):
			base = self.consensus_at_pos(seq_name, position)
			
			consensus += base #add base 
			
			if base not in ["", "N"]:
				mapped = True
		
		#2 no reads map return clause
		if mapped == False:
			return ''
		else:
			return consensus

	def best_consensus(self):
		"""
		return the majority base in reference with best mapping
		best mapping - more positions = better mapping
		"""
		references = []
		best_sequence = 'x'
		best_map = -1

		#obtain all mapped reads that do not repeat
		for read in self.reads:
			if read.is_mapped and read.is_primary:
				if read.rname not in references:
					references.append(read.rname)
	
		for ref in references: 
			cov = 0
			ref_len = self.references.get(ref, 0)

			for i in range(1, ref_len + 1):
				check = [c for c in self.reads if c.is_mapped and c.is_primary and c.rname == ref and c.pos<= i < c.pos + c.mapped_len and c.read_idx_at_pos(i) != []]

				if len(check) > 0:
					cov += 1

			if cov > best_map: #update for best 
				best_map = cov
				best_sequence = ref
		
		if best_sequence == "" or best_map < 0:
			return ""
		
		return(self.consensus(best_sequence))

