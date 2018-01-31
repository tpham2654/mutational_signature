import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def reversecomplement(dnaseqin):
	dnain = Seq(dnaseqin,generic_dna)
	return str(dnain.reverse_complement())

def makedictkey(change,sthreebefore):
	arrchange = change.split(">")
	nucbefore = arrchange[0]
	nucafter = arrchange[1]
	threebefore = [s for s in sthreebefore]
	threeafter = [s for s in threebefore]
	if threebefore[1] == nucbefore:
		threeafter[1] = nucafter
		return str("".join(threebefore)) + " ==> " + str("".join(threeafter))
	return False

def makereversecomplementdictkey(dictkeyin):
	arrdictkey = dictkeyin.split(" ==> ")
	threebefore = reversecomplement(arrdictkey[0].strip())
	threeafter = reversecomplement(arrdictkey[1].strip())
	return threebefore + " ==> " + threeafter
	

def makesigdict(signum):
	mutdf = pd.read_table('signatures_probabilities.txt')
	dictout = {}
	for i,r in mutdf.iterrows():
		dictkey = makedictkey(r["Substitution Type"],r["Trinucleotide"])
		revcomplementdictkey = makereversecomplementdictkey(dictkey)
		if dictkey != False:
			dictout[dictkey] = r["Signature " + str(signum)]
			dictout[revcomplementdictkey] = r["Signature " + str(signum)]
	return dictout
