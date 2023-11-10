import sys
from getOptPar import getOptPar
import numpy as np


def eliminateRedundancy(fasta_db, new_fasta_db, redundant_ids_db):

	fasta_db = [x.rstrip() for x in open(fasta_db).readlines()]
	seq_to_protid_dict = {}
	prot_ids = []
	seqs = []
	for line in fasta_db:
		if line[0] == '>':
			prot_id = line[1::]
			prot_ids.append(prot_id)
		else:
			seq = line
			seqs.append(seq)
			if seq not in seq_to_protid_dict:
				seq_to_protid_dict[seq] = []
			seq_to_protid_dict[seq].append(prot_id)


	uniq_seqs = list(set(seqs))

	nr_fasta_file = open(new_fasta_db, 'w')
	redundancy_file = open(redundant_ids_db, 'w')
	for seq in uniq_seqs:
		redundancy_file.write(seq_to_protid_dict[seq][0].replace(' ', '')+'\t'+','.join(seq_to_protid_dict[seq]).replace(' ', '')+'\n')
		nr_fasta_file.write('>'+seq_to_protid_dict[seq][0].replace(' ', '')+'\n')
		nr_fasta_file.write(seq+'\n')
		



	


def eliminateRedundancyWindow(fasta_db, w, new_fasta_db):



	fasta_db = [x.rstrip() for x in open(fasta_db).readlines()]
	db_dict = {}
	labels = []
	for line in fasta_db:
		if line[0] == '>':
			db_dict[line] = []
			key = line
			labels.append(key)
		else:
			db_dict[key].append(line)
			db_dict[key].append(len(line))
			db_dict[key].append(splitString(line, w))
		


	# sort descendently by length
	lens = np.array([db_dict[x][1] for x in labels])
	labels = np.array(labels)

	o = lens.argsort()[::-1]
	lens = lens[o]
	labels = labels[o]



	stop = False
	maintain = np.ones((len(labels),))
	for i in range(len(labels)):
		if maintain[i]:
			windowsi = db_dict[labels[i]][2]
			for j in range(i+1, len(labels)):
				if maintain[j]:
					windowsj = db_dict[labels[j]][2]
					overlap = checkOverlap(windowsi, windowsj)
					if overlap:
						maintain[j] = 0
	labels = labels[np.array(maintain, dtype=bool)]


	f = open(new_fasta_db, 'w')	
	for li in labels:
		f.write(li+'\n')
		f.write(db_dict[li][0]+'\n')

	return db_dict	




def splitString(s, w):
	l = []
	n = len(s)
	for i in range(n-w):
		l.append(s[i:i+w+1])	
	return l

def checkOverlap(windowsi, windowsj):
	overlap = False
	for wi in windowsi:
		for wj in windowsj:
			if wi == wj:
				overlap = True
				return overlap
	return overlap
		



if __name__ == '__main__':
	iargs = sys.argv
	fasta_db = iargs[1]
	w = getOptPar(iargs, 'w', 9)
	eliminateRedundancy(fasta_db, w)
