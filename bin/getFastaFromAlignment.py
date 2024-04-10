import sys
from getOptPar import getOptPar
from fasta_utilities import read_fasta

def filter_alignment(align_file, outfile, fasta_db, filter_par, filter_value, W):

	fasta_db = read_fasta(fasta_db)
	fasta_db_new = {}
	for key in fasta_db:
		new_key = key.split(' ')[0].replace('>', '')
		fasta_db_new[new_key] = fasta_db[key]
	fasta_db = fasta_db_new

	align_file = [x.rstrip().split(',') for x in open(align_file).readlines()]
	out_cont = []
	for line in align_file:
		subjct = line[4] #human
		subject_aligned_seq = line[5]
		subject_start = line[6]
		subject_end = line[7]
		query_aligned_seq = line[1]

		if '-' in subject_aligned_seq:
			continue
		start = int(line[6])-1
		end = int(line[7])


		query = '|'.join([line[0], line[2], line[3]]) #infct
		ident = str(round(float(line[9]),1))
		evalue = line[13]
		bitscore = line[8]


		
		acceptable_score = False
		if filter_par == 'bitscore':
			score = float(line[8])
			if score > filter_value:
				acceptable_score = True
			else:
				acceptable_score = False

		elif filter_par == 'evalue':
			score = float(line[13])
			if score < filter_value: 
				acceptable_score = True
			else:
				acceptable_score = False
		
		elif filter_par == 'identity':
			score = float(line[9])
			if score > filter_value: 
				acceptable_score = True
			else:
				acceptable_score = False
		
		

		if not acceptable_score:
			continue
		
		
		whole_seq = fasta_db[subjct]
		aligned_seq = whole_seq[start:end]
		#print([len(aligned_seq), aligned_seq])

		len_seq = len(aligned_seq)


		if len_seq >= W:
			output_seq = aligned_seq
		else:
			diff = W-len_seq
			if diff%2==0: #len_seq is pari 
				start_  = start - int(diff/2)				
				end_  = end + int(diff/2)
					
			else: #len_seq is dispari (shift to left side)
				start_  = start - int(diff/2+0.5)				
				end_  = end + int(diff/2-0.5)	
				
			if start_ < 0:
				end_ = end_ + (-start_)
				start_ = 1
				Wmer_seq = whole_seq[0:end_]	

			elif end_ > len(whole_seq)-1:
				start_ = start_ - (end_-len(whole_seq)) 
				end_ = len(whole_seq)				
				Wmer_seq = whole_seq[start_:None]
				
			else:
				Wmer_seq = whole_seq[start_:end_]

			start = start_
			end = end_
			output_seq = Wmer_seq

		
			"""
			else:

				out_cont.extend(['>'+epitope_name, whole_seq])
				print(['>'+epitope_name, whole_seq])


				#print(len(Wmer_seq))		
				out_cont.extend(['>'+epitope_name, Wmer_seq])
			"""

		epitope_name = '|'.join([subjct, str(start+1), str(end)])
		epitope_name += '@@@'+query
		epitope_name += '@@@'+'|'.join([ident, evalue, bitscore, subject_aligned_seq, subject_start, subject_end, query_aligned_seq])+'@@@'
		out_cont.extend(['>'+epitope_name, output_seq])

	if len(out_cont) == 0:
		return False
	else:
		f = open(outfile, 'w')
		f.write("\n".join(out_cont)+"\n")
		return True




	

if __name__ == '__main__':
	iargs = sys.argv
	align_file = iargs[1]
	outfile = iargs[2]
	fasta_db = iargs[3]
	filter_par = getOptPar(iargs, 'filter_par', 'evalue', 'str')
	filter_value = getOptPar(iargs, 'filter_value', 0.05, 'float')
	W = getOptPar(iargs, 'W', 15, 'int')
	filter_alignment(align_file, outfile, fasta_db, filter_par, filter_value, W)
