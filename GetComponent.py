#!/etc/bin/env python
# -*- coding: utf-8 -*-
#

import os
import pickle
import numpy as np
from subprocess import call

AminoAcid = ['G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','D','E','K','R','H']
with open('./Dataset/Dipeptide', 'r') as f:
	Dipeptide = pickle.load(f)

with open('./Dataset/Tripeptide', 'r') as f:
	Tripeptide = pickle.load(f)

def checkSeq(sequence):
	pass

def AminoAcid2Freq(sequence):
	""" 
	Calculate the frequency of 20 amino acids in one Sequence (string)
	type sequence: string
	rtype: list[float]
	"""
	seq2upper = sequence.upper()
	freq = []
	for i in range(20):
		freq.append(seq2upper.count(AminoAcid[i])/float(len(seq2upper)))
	return (AminoAcid, freq)

def Dipeptide2Freq(sequence, gap=0):
	"""
	Calculate the frequency of 400 depeptides in one Sequence (string)
	type sequence: string
	type gap: int
	rtype: list[float]
	"""
	seq2upper = sequence.upper()
	counts = [0]*400
	for i in range(len(seq2upper)-gap-1):
		for m in range(20):
			if seq2upper[i] == AminoAcid[m]:
				break
		for n in range(20):
			if seq2upper[i+gap+1] == AminoAcid[n]:
				break
		counts[20*m+n] += 1
	freq = [0]*400
	for i in range(400):
		freq[i] = counts[i]/float(len(seq2upper)-gap-1)
	return (Dipeptide, freq)

def Tripeptide2Freq(sequence):
	"""
	Calculate the frequency of 8000 Tripeptides in one Sequence (string)
	type sequence: string
	rtype: list[float]
	"""
	seq2upper = sequence.upper()
	counts = [0]*8000
	for i in range(len(seq2upper)-2):
		for l in range(20):
			if seq2upper[i] == AminoAcid[l]:
				break		
		for m in range(20):
			if seq2upper[i+1] == AminoAcid[m]:
				break
		for n in range(20):
			if seq2upper[i+2] == AminoAcid[n]:
				break
		counts[400*l+20*m+n] += 1 
	freq = [0]*8000
	for i in range(8000):
		freq[i] = counts[i]/float(len(seq2upper)-2)
	return (Tripeptide, freq)

def svmfile(obj, filename, *args):
	'''
	Output file with libsvm data type
	'''
	batch = 0
	with open(filename, 'w') as f:
		for lab, arg in zip(range(len(args)), args):
			for i in range(arg):
				f.write('{} '.format(lab))

				idx = i + batch
				for j, ele in zip(range(1, 8001), obj[idx]):
					f.write('{}:{} '.format(j, ele))

				f.write('\n')

			batch += arg

def gridfile(filename):
	'''
	Read data.grid, and return freq
	'''
	freq = {}
	with open(filename, 'r') as f:
		lines = [line.strip() for line in f]

	items = lines[-1].split(' ')
	for idx, item in zip(('c', 'g', 'rate'), items):
		freq[idx] = float(item)

	return freq

def svmfunc():
	pass

if __name__ == '__main__':
	
	m = 178
	n = 226
	os.chdir("/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_CancerLectin/")
		
	Lectin = []
	# process CancerLectin and nonCancerLectin
	with open("./Dataset/cancerlectin178.txt", 'r') as f:
		seqs = [line.strip() for line in f if line[0] != '>']

	for seq in seqs:
		t, v = Tripeptide2Freq(seq)
		Lectin.append(v)
		
	with open("./Dataset/non-cancerlectin226.txt", 'r') as f:
		seqs = [line.strip() for line in f if line[0] != '>']
		
	for seq in seqs:
		t, v = Tripeptide2Freq(seq)
		Lectin.append(v)
	
	Lectin = np.array(Lectin, dtype=float)
	#svmfile(Lectin, './Dataset/data.txt', m, n)
	'''
	with open('./Dataset/Lectin', 'w') as f:
		pickle.dump(Lectin, f)
	
	with open('./Dataset/Lectin', 'r') as f:
		Lectin = pickle.load(f)
	'''

	Book = {}
	MF = np.empty(shape=(404, 1), dtype=float)
	orderBook = []
	while len(Book) < 8000:
		Freq = {}
		for idx in range(8000):
			if np.sum(Lectin[:, idx]) < 0.0001:
				'0.0001: 1 / 10000'
				Book[idx] = False
				continue

			if not idx in Book:
				A = np.insert(MF[:, :-1], 0, Lectin[:, idx], axis=1)
				# call libsvm
				svmfile(A, './tmp/data.ori', m, n)
				call(('svm-scale -l 0 -u 1 ./tmp/data.ori > ./tmp/data.scale'), shell=True)
				call(('grid.py -v 4 ./tmp/data.scale > ./tmp/data.grid'), shell=True)
				freq = gridfile('./tmp/data.grid')
				
				Freq[idx] = freq

		mrate = {'c':0.0, 'g':0.0, 'rate':0.0}
		mkey = 0
		for key, val in Freq.items():
			if val['rate'] > mrate['rate']:
				mrate = val
				mkey = key

		orderBook.append(mkey)
		Book[mkey] = (len(orderBook), mrate)
		print 'Processing [{}] features, the predictive rate is: [{}]'.format(len(orderBook), mrate['rate'])

		# Save the results
		with open('result', 'w') as f:
			pickle.dump((orderBook, Book), f)
	
		MF = np.insert(MF, 0, Lectin[:, mkey], axis=1)
			
	
