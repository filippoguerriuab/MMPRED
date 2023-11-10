import sys
from getOptPar import getOptPar

def read_fasta(fasta, dict_or_lists = 'dict'):

	f = [line.rstrip() for line in open(fasta).readlines()]
	fasta_dict = {}
	for line in f:
		if line[0] == '>':
			id = line
			fasta_dict[id] = ''
		else:
			fasta_dict[id] += line


	if dict_or_lists == 'dict':
		return fasta_dict

	elif dict_or_lists == 'lists':
		seq_name = []
		seq_seq = []
		for id in fasta_dict:
			seq_name.append(id)
			seq_seq.append(fasta_dict[id])
		return seq_name, seq_seq


def filter_by_length(fasta, l, sign):



	if type(fasta) == str:

		fasta_dict = read_fasta(fasta)



	if sign == 'gte':
		fasta_dict_filtered = {key : fasta_dict[key] for key in fasta_dict if len(fasta_dict[key]) >= l}
	elif sign == 'lwe':
		fasta_dict_filtered = {key : fasta_dict[key] for key in fasta_dict if len(fasta_dict[key]) <= l}
	
	for ids in fasta_dict_filtered:
		print(ids)
		print(fasta_dict[ids]) 

def filter_by_proteinset(fasta, prot_set, include_exclude='include'):

	fasta_dict = read_fasta(fasta, 'dict')
	prot_set = set([x.rstrip() for x in open(prot_set).readlines()])




	ids_matching = []
	for key in fasa_dict:
		uniprot_id = key.split('|')[0].replace('>', '')
		if uniprot_id in prot_set:
			ids_matching.append(key)



	if include_exclude == 'exclude':
		#print(['a', len(fasta_dict)])
		for idi in ids_matching:
			del fasta_dict[idi]
		#print(['b', len(fasta_dict)])


	elif include_exclude == 'include':
		fasta_dict = {idi:fasta_dict[idi] for idi in ids_matching}


	for ids in fasta_dict:
		print(ids)
		print(fasta_dict[ids])




def filter_by_ids(fasta, ids_list, include_exclude='include'):

	fasta_dict = read_fasta(fasta, 'dict')
	ids_list = set([x.rstrip() for x in open(ids_list).readlines()])




	ids_matching = []
	for key in fasta_dict:
		if key in ids_list:
			ids_matching.append(key)



	if include_exclude == 'exclude':
		#print(['a', len(fasta_dict)])
		for idi in ids_matching:
			del fasta_dict[idi]
		#print(['b', len(fasta_dict)])


	elif include_exclude == 'include':
		fasta_dict = {idi:fasta_dict[idi] for idi in ids_matching}


	for ids in fasta_dict:
		print(ids)
		print_prot_seq(fasta_dict[ids])


def split_by_w(fasta, W):
	fasta_dict = read_fasta(fasta)
	for seq_id in fasta_dict:
		seq = fasta_dict[seq_id]
		for i in range(len(seq)-W+1):
			curr_seq_id = '|'.join([seq_id, str(i+1), str(i+W)])
			if 'X' not in seq[i:i+W]:
				print(curr_seq_id)
				print(seq[i:i+W])


	



def print_prot_seq(seq):
	W = 80
	[print(seq[i:i+W]) for i in range(0,len(seq),W)]






if __name__ == '__main__':
	fasta = sys.argv[1]
	option = sys.argv[2]
	if option == 'filter_by_length':
		sign = sys.argv[3]
		l = int(sys.argv[4])
		filter_by_length(fasta, l, sign)
	elif option == 'filter_by_proteinset':
		include_exclude = sys.argv[3]
		prot_set = sys.argv[4]
		filter_by_set(fasta, prot_set, include_exclude)
	elif option == 'filter_by_ids':
		include_exclude = sys.argv[3]
		ids_list = sys.argv[4]
		filter_by_ids(fasta, ids_list, include_exclude)
	elif option == 'split_by_W':
		W = int(sys.argv[3])
		split_by_w(fasta, W)
	