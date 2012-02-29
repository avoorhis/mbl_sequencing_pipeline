#!/usr/bin/python
import string
import sys
import os
import copy
import types

"""
Author: Pierre Tuffery, 2004-2008
Version 3.0
Ressource Parisienne en Bioinformatique Structurale
http://bioserv.rpbs.univ-paris-diderot.fr

This is free software. You can use it, modify it, distribute it.
However, thanks for the feedback for any improvement you  bring to it!

Simple classes to manage fasta, multifasta data.

Example:
from Fastav3 import *
x = fasta("toto.fst")
y = fasta("toto2.fst")
z = x+y
k = z - y
k == x

a = x['d1dlwa_']
b = x['d1dlya_']
a in a
c = a+b
d = c-b
d == a

"""

##
## Class related to one single sequence
##
class sfasta:
	"""
	sfasta: a class to manage a single fasta sequence
	data is organized as a dictionary of:
	id: the sequence identifier
	cmt: (comment after id)
	s: the sequence
	"""
	def __init__(self, id = None, seq = None, cmt = None, verbose = 0):
		"""
		id: sequence id
		cmt: comment on the "> id" line, after the id
		seq: the sequence itself
		"""
		self.data = {"s":seq, "cmt": cmt, "id": id}

	def __getitem__(self,pos):
		return self.data["s"][pos]

	def __repr__(self):
		"""
		Flat representation of sequence
		"""
		rs = ""
		rs += ">%s %s\n" % (self.data["id"], self.data["cmt"])
		rs += "%s\n" % self.data["s"]
		return rs

	def __len__(self):
		"""
		return sequence length
		"""
		return len(self.data["s"])

	def __contains__(self, other):
		"""
		sfasta.__contains__(other) :
		does sequence contain some subsequence ?
		if x in y: # is a a subsequence of y ?
		"""
		if other.data["s"] in self.data["s"]:
			return True
		return False

	def __add__(self, other):
		"""
		sfasta.__add__ (other)
		concatenate two sequences in a new one.
		c = a + b # c is sequence of a then sequence of b merged into one
		"""
		out = sfasta()
		out.data["s"] = self.data["s"] + other.data["s"]
		out.data["id"] = "%s+%s" % (self.data["id"],other.data["id"])
		out.data["cmt"] = "%s+%s" % (self.data["cmt"], other.data["cmt"])
		return out

	def __sub__(self, other):
		"""
		sfasta.__sub__(other) :
		remove exact occurrence of other in sequence ?
		c = a - b # c is sequence of a from which b has been removed
		"""
		out = sfasta()
		out.data["s"] = self.data["s"]
		out.data["id"] = self.data["id"]
		out.data["cmt"] = self.data["cmt"]
		if other in self:
			pos = self.data["s"].index(other.data["s"])
			out.data["s"] = self.data["s"][:pos] + self.data["s"][pos+len(other.data["s"]):]
			if ("-%s" % other.data["id"]) in out.data["id"]:
				seed = "-%s" % other.data["id"]
			elif ("+%s" % other.data["id"]) in out.data["id"]:
				seed = "+%s" % other.data["id"]
			elif ("%s-" % other.data["id"]) in out.data["id"]:
				seed = "%s-" % other.data["id"]
			elif ("%s+" % other.data["id"]) in out.data["id"]:
				seed = "%s+" % other.data["id"]
			else:
				seed = None
			if seed:
				pos = self.data["id"].index(seed)
				out.data["id"] = self.data["id"][:pos] + self.data["id"][pos+len(seed):]
			else:
				out.data["id"] = "%s-%s" % (self.data["id"],other.data["id"])


			if ("-%s" % other.data["cmt"]) in out.data["cmt"]:
				seed = "-%s" % other.data["cmt"]
			elif ("+%s" % other.data["cmt"]) in out.data["cmt"]:
				seed = "+%s" % other.data["cmt"]
			elif ("%s-" % other.data["cmt"]) in out.data["cmt"]:
				seed = "%s-" % other.data["cmt"]
			elif ("%s+" % other.data["cmt"]) in out.data["cmt"]:
				seed = "%s+" % other.data["cmt"]
			else:
				seed = None


			if seed:
				pos = self.data["cmt"].index(seed)
				out.data["cmt"] = self.data["cmt"][:pos] + self.data["cmt"][pos+len(seed):]
			else:
				out.data["cmt"] = "%s-%s" % (self.data["cmt"],other.data["cmt"])
		return out

	def __eq__(self, other):
		"""
		sfasta.__eq__() :
		are the sequences identical ?
		a == b # are sequences strictly identical ?
		"""
# 		if len(self) != len(other):
# 			return False
		if other == None:
			return False
		if self.data["s"] == None:
			return False
		if self.data["s"] != other.data["s"]:
			return False
		return True

	def __ne__(self, other):
		"""
		sfasta.__ne__() :
		are the sequences identical ?
		a != b # are sequences not strictly identical ?
		"""
		if len(self) != len(other):
			return True
		if self.data["s"] != other.data["s"]:
			return True
		return False

	def s(self, seq = None):
		"""
		return the sequence string (if seq is None), or assign it (if seq is specified)
		"""
		if seq == None:
			return self.data["s"]
		else:
			self.data["s"] = seq

	def cmt(self, cmt = None):
		"""
		return comment, or assign it
		"""
		if cmt == None:
			return self.data["cmt"]
		else:
			self.data["cmt"] = cmt

	def id(self, id = None):
		"""
		return id, or assign it
		"""
		if id == None:
			return self.data["id"]
		else:
			self.data["id"] = id

	def out(self, f = sys.stdout, oneLine = True, step = 80, upper = False, lower = False, star = False, pretty = False):
		"""
		A sequence formatter.

		out: will output formatted content of sequence
		oneLine: all sequence on one line, else break at step
		step: if oneLine is False, line are truncated each step.
		upper: force uppercase
		lower: force lowercase
		star: add star at end of sequence
		pretty: split each line as series of 10 letters separated by blank
		"""
		if self.data["cmt"] == None or self.data["cmt"] == '':
		    rs = ">%s" % (self.data["id"],)
		else:
		    rs = ">%s %s" % (self.data["id"], self.data["cmt"])
		if oneLine:
			step = len(self.data["s"]) + 10
		seqData = self.data["s"]
		outlist = []
		while len(seqData) > 0:
			if upper:
				outlist.append("%s" % seqData[:step].upper() )
			elif lower:
				outlist.append("%s" % seqData[:step].lower() )
			else:
				outlist.append("%s" % seqData[:step] )
			seqData = seqData[step:]
		for l in outlist:
			ol = l
			if pretty:
				ol = ""
				while len(l) > 0:
					ol += "%s " % l[:10]
					l = l[10:]
			rs += "\n%s" % ol
		if star:
			rs += "*"
		rs += "\n"
		f.write(rs)
			
	def write(self,fname, fmode = "w", oneLine = True, step = 80, upper = False, lower = False, star = False, pretty = False):
		"""
		write: this will perform sequence output in fname, using fmode (one of classical "w", "a", etc)
		"""
		f = open(fname,fmode)
		self.out(f, oneLine, step, upper, lower, star)
		f.close()

##
## Class related to fasta formated sequence file(s)
##
class fasta:
	"""
	a class to multi fasta data
	each sequence is a instance of the sfasta class
	
	Various operators are defined to merge, combine fasta collections and subsets.
	"""

	def __init__(self, fname = "", setName = None, verbose = 0):
		"""
		fname: if fname is a string, then it is assimilated to a filename
		if it is a list, then it is assimilated to a list of lines to parse.
		setName: a name for the set.

		data = fasta("toto.fst")
		"""
		
		if isinstance(fname,types.StringType):  # read file from disk
			self.data = {}
			self.nSeq = 0
			self.name = fname
			if fname != "":
				self.load(fname,verbose)
		elif isinstance(fname,types.ListType):    # already a list of lines
			if setName != None:
				self.name = setName
			else:
				self.name = "stdin"
			self.data = {}
			self.nSeq = 0
			self.parse(fname, verbose)
		
	def load(self,fname, verbose = 0):
		"To load a generic fasta file from disk."
		lines = open(fname,'r').readlines()
		if verbose:
			print fname, ": Read ",len(lines),' lines'
		self.parse(lines, verbose)

	def parse(self, lines, verbose):
		""" 
		perform the effective parsing of lines
		(i.e. a series of lines as:
		> Id comment OR >Id comment
		dataline
		dataline

		Each sequence is a dictionnary of
		id
		comment
		sequence
		"""

		# Split all lines so a to remove blanks
		# for i in range(0,len(lines)):
		#	lines[i] = string.split(lines[i])

		self.data = {}
		self.nSeq = 0
		curId = None
		for l in lines:
			it = l.split()
			if it[0][0] == ">":

				# Add previous data
				if curId != None:
					# self.data[curId] = {"s":seq, "cmt": comment, "id": curId}
					self.data[curId] = sfasta(id = curId, seq = seq, cmt = comment)
				curId = None
				# parse current data
				self.nSeq += 1
				# Find the Id
				if it[0] == ">":
					curId = it[1]
					try:
						comment = " ".join(it[2:])
					except:
						comment = ""
				else:
					curId = it[0][1:]
					try:
						comment = " ".join(it[1:])
					except:
						comment = ""
				seq = ""
			else:
				seq = seq + "".join(it)

		# Add previous data
		if curId != None:
			# self.data[curId] = {"s":seq, "cmt": comment, "id": curId}
			self.data[curId] = sfasta(id = curId, seq = seq, cmt = comment)
			
		if verbose:
			print "%s : Read %d sequences" % (self.name,self.nSeq)

	def __getitem__(self,key):
		"""
		accessor to one sequence: data["myId"]
		"""
		return self.data[key]

	def __delitem__(self,key):
		"""
		destruction of one sequence: del data["myId"]
		"""
		del self.data[key]
		self.nSeq -= 1

	def __setitem__(self, id , seq):
		"""
		insertion of one sequence: data["myId"] = seq # where seq is a sfasta instance
		"""
		self.data[id] = seq
		self.nSeq += 1
	
	def __contains__(self, other):
		"""
		 x in y # where the match is based on sequence Ids
		"""
		oIds = other.ids()
		sIds = self.ids()
		for aId in oIds:
			if aId not in sIds:
				return False
		return True

	def __add__(self, other):
		"""
		fasta.__add__ (other)
		concatenate two multifasta objects in a new one
		"""
		outSet = fasta()

##		outSet.data = self.data[:]
		for aId in self.ids():
			outSet[aId] =  self.data[aId]
		for aId in other.ids():
			outSet[aId] =  other.data[aId]
		return outSet

	def __sub__(self, seq):
		"""
		 c = x - y # remove sequences of y present in x, return in c
		"""
		outSet = fasta()
		for aId in self.ids():
			outSet[aId] =  self.data[aId]
		for aId in seq.ids():
			outSet.__delitem__(aId)
		return outSet

	def __or__(self,other): # self | other : union
		"""
		merge two multi fasta objects, based on sequence Ids
		"""
		outSet = fasta()

##		outSet.data = self.data[:]
		for aId in self.ids():
			outSet[aId] =  self.data[aId]
		for aId in other.ids():
			if aId not in outSet.ids():
				outSet[aId] = other.data[aId]
		return outSet
	
	def __and__(self,other):  # self & other : intersection
		"""
		intersection of two multifasta (based on Ids)
		"""
		outSet = fasta()

		for aId in self.ids():
			if aId in other.ids():
				outSet[aId] = self.data[aId]
		return outSet

	def __eq__(self, other):
		"""
		x == y # True is x and y correspond to the same collections of Ids
		e.g.:
		c = a + b
		d = c - a
		d == b # True
		"""
		if len(self) != len(other):
			return False
		sIds = self.ids()
		for aId in other.ids():
			if aId not in sIds:
				return False
		return True
	
# 	def __add__(self,other):
# 		return self | other
	def __len__(self):
		"""
		The number of sequences if the collection
		"""
		return len(self.ids())

	def __repr__(self):  # pour print
		res = ""
		for aId in self.ids():
			res = res + self.data[aId].__repr__()
		return res				

	def ids(self):
		"""
		fasta.ids:
		return a list of Ids of the sequences
		"""
		return self.data.keys()

##
## write sequence(s) to stdout
##
		
	def out(self, f = sys.stdout, Ids = None, oneLine = True, step = 80, upper = False, lower = False, star = False, pretty = False):
		"""
		A sequence formatter.

		f: the file descriptor to write the data content
		Ids: a selection of Ids to write. If None: everything is output.
		oneLine: all sequence on one line, else break at step
		step: if oneLine is False, line are truncated each step.
		upper: force uppercase
		lower: force lowercase
		star: add star at end of sequence
		pretty: split each line as series of 10 letters separated by blank
		"""

		if Ids == None:
			Ids = self.ids()
		for aId in Ids:
			self.data[aId].out(f = f, oneLine = oneLine, step = step, upper = upper, lower = lower, star = star, pretty = pretty)

	def write(self,fname, fmode = "w", Ids = None, oneLine = True, step = 80, upper = False, lower = False, star = False, pretty = False):
		"""
		fasta.write:
		write a collection (subset or all) of sequences to file
		fname: filename
		fmode: one of "w", "a", etc
		Ids: if None all sequences are output. Else, only the Ids in the list are output
		will propagate attribute parameters to fasta.out()
		"""
		f = open(fname, fmode)
		self.out(f = f, Ids = Ids, oneLine = oneLine, step = step, upper = upper, lower = lower, star = star, pretty = pretty)
		f.close()
				

## #
## # A composite fasta file into series of simple fasta files
## #
	def splitwrite(self,fileExt=".fst", path="./", Ids = None, oneLine = True, step = 80, upper = False, lower = False, star = False, pretty = False):
		"""
		fasta.splitwrite():
		This will split output on the form one file per sequence.
		path: the directory to write in (./)
		fileExt: file extension to use (.fst)
		Ids: if not None, only these sequences will be output.
		"""
		if path[-1] != "/":
			path = path+"/"

		if Ids == None:
			Ids = self.ids()
		if path[-1] != '/':
			path = path+"/"
		for aId in Ids:
			self[aId].write("%s%s%s" % (path, aId, fileExt),oneLine = oneLine, step = step, upper = upper, lower = lower, star = star, pretty = pretty)


	def subSet(self, theList, verbose = 0):
		"""
		return a new instance corresponding to the subset of Ids in theList.
		"""
		outSet = fasta()

		for aId in theList:
			print aId
			outSet.data[aId] = self.data[aId]
		if verbose:
			print "Sub selection of ",len(outSet.data),"(wanted ",len(theList),") among ",len(self.data)
		return outSet


if __name__=='__main__':
	x = fasta("toto.fst")
	x
