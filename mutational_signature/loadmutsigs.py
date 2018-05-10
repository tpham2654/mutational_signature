import pandas as pd


def reversecomplement(dna_seq_in):
	dna_complement = ""
	complement_pair = {"A":"T",
	"G":"C",
	"T":"A",
	"C":"G"}
	for nt in dna_seq_in:
		dna_complement = dna_complement + complement_pair[nt]
	return str(dna_complement[::-1]) #reversed

def makedictkey(change,sthreebefore):
	#changes COSMIC formatted mutation description to dictionary key format this program uses.
	arr_nucleotidechange = change.split(">")
	nucleotide_before = arr_nucleotidechange[0]
	nucleotide_after = arr_nucleotidechange[1]
	arr_three_nt_before = [s for s in sthreebefore]
	arr_three_nt_after = [s for s in arr_three_nt_before]
	if arr_three_nt_before[1] == nucleotide_before:
		arr_three_nt_after[1] = nucleotide_after
		return str("".join(arr_three_nt_before)) + " ==> " + str("".join(arr_three_nt_after))
	return False

def makereversecomplementdictkey(dictkeyin):
	arr_dictkey = dictkeyin.split(" ==> ")
	three_nt_before = reversecomplement(arr_dictkey[0].strip())
	three_nt_after = reversecomplement(arr_dictkey[1].strip())
	return three_nt_before + " ==> " + three_nt_after
	

def makesigdict(signum,doreversecomplement=True):
	mutation_sigs_df = pd.read_table('signatures_probabilities.txt') #file comes from COSMIC.
	mutsig_dictout = {}
	for i,r in mutation_sigs_df.iterrows():
		mutsig_dictkey = makedictkey(r["Substitution Type"],r["Trinucleotide"])
		revcomplement_dictkey = makereversecomplementdictkey(mutsig_dictkey)
		if mutsig_dictkey != False:
			mutsig_dictout[mutsig_dictkey] = r["Signature " + str(signum)]
			#reverse complement of described mutation is treated the same as original.
			if doreversecomplement==True:
				mutsig_dictout[revcomplement_dictkey] = r["Signature " + str(signum)] 
	return mutsig_dictout
	
def subtract_sigs(sig_to_subtract_from_num,sig_to_subtract_with_num):
	sig_out = dict()
	sig_to_subtract_from = makesigdict(sig_to_subtract_from_num)
	sig_to_subtract_with = makesigdict(sig_to_subtract_with_num)
	sig1_in_2 = [k for k in sig_to_subtract_from if k in sig_to_subtract_with]
	sig2_in_1 = [k for k in sig_to_subtract_with if k in sig_to_subtract_from]
	common_sigs = list(set([k for k in sig1_in_2+sig2_in_1]))
	only_in_sig_to_subtract_from = [k for k in sig_to_subtract_from if k not in sig_to_subtract_with]
	only_in_sig_to_subtract_with = [k for k in sig_to_subtract_with if k not in sig_to_subtract_from]
	for key in common_sigs:
		subtracted_probability = sig_to_subtract_from[key] - sig_to_subtract_with[key]
		if(float(subtracted_probability)>float(0)):
			sig_out[key] = subtracted_probability
		else:
			sig_out[key] = 0
	for key in only_in_sig_to_subtract_from:
		sig_out[key] = only_in_sig_to_subtract_from[key]
	for key in only_in_sig_to_subtract_with:
		sig_out[key] = only_in_sig_to_subtract_with[key]
	return sig_out