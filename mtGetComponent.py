#!/etc/bin/env python
# -*- coding: utf-8 -*-
#

import os
import pickle
import numpy as np
from atexit import register
from threading import Thread, Lock, current_thread, active_count
from time import ctime
from subprocess import call
from copy import deepcopy

AminoAcid = ['G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','D','E','K','R','H']
with open('./Dataset/Dipeptide', 'r') as f:
	Dipeptide = pickle.load(f)

with open('./Dataset/Tripeptide', 'r') as f:
	Tripeptide = pickle.load(f)

lock = Lock()
lockFreq = Lock()
Freq = {}

class cleanOutputSet(set):
	def __str__(self):
		return ', '.join(x for x in self)

remaining = cleanOutputSet()

class MyThread(Thread):
	def __init__(self, func, obj, fold, idx, name='mythread', *args):
		Thread.__init__(self)
		self.__func = func
		self.__obj = obj
		self.__fold = fold
		self.__idx = idx
		self.__args = args
		self.__name = name

	def run(self):
		self.__func(self.__obj, self.__fold, self.__idx, self.__name, self.__args)

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

def checkSeq(sequence):
	pass

def svmfunc(obj, fold, idx, myname, args):
	with lock:
		remaining.add(myname)
		print '[{}] Started {}'.format(ctime(), myname)

	svmfile(obj, './tmp/{0}.ori'.format(myname), args)
	call(('svm-scale -l 0 -u 1 ./tmp/{0}.ori > ./tmp/{0}.scale'.format(myname)), shell=True)
	call(('grid.py -v {0} ./tmp/{1}.scale > ./tmp/{1}.grid'.format(fold, myname)), shell=True)
	freq = gridfile('./tmp/{0}.grid'.format(myname))

	with lock:
		remaining.remove(myname)
		print '[{}] Comleted {}'.format(ctime(), myname)
		print '    {}'.format(freq)
		print '    (remaining: {})'.format(remaining or None)
	
	with lockFreq:
		Freq[idx] = freq

def svmfile(obj, filename, args):
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

def rmNoSenceFeature(data, Book):
	nCol = data.shape[1]
	for idx in range(nCol):
		try:
			if (not idx in Book) and (np.sum(data[:, idx] ) < 0.0001):
				'0.0001: 1 / 10000'
				Book[idx] = False
				
		except ValueError:
			pass

	return Book

#if __name__ == '__main__':
def optimalFunc(Lectin, Book={}, orderBook=[], MF=np.empty(shape=(404, 1), dtype=float), nthreads = 2):
		
	m = 178
	n = 226
	nthreads = nthreads
	nCol = Lectin.shape[1]
	Book = rmNoSenceFeature(Lectin, Book)
	while len(Book) < nCol:
		Freq.clear()
		flag = 0
		threads = []

		NUM = nCol - len(Book)
		nCir = NUM / nthreads
		cir = 1
		for idx in range(nCol):
			if not idx in Book:
				A = np.insert(MF[:, :-1], 0, Lectin[:, idx], axis=1)
				# call libsvm
				threads.append(MyThread(svmfunc, A, 4, idx, 'Thread-{}'.format(flag+1), m, n ) )
				flag += 1
				
				if (cir <= nCir) and (flag == nthreads):
					for i in range(flag):
						threads[i].start()
					for i in range(flag):
						threads[i].join()

					#print '{} thread(s) is alive'.format(active_count())
					cir += 1
					flag = 0
					threads = []
					
				if (cir == nCir + 1) and (flag == NUM % nthreads):
					for i in range(flag):
						threads[i].start()
					for i in range(flag):
						threads[i].join()				

					print '{} thread(s) is alive'.format(active_count())
					cir += 1
					flag = 0
					threads = []
		
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
		with open(r'result', 'w') as f:
			pickle.dump((orderBook, Book), f)
	
		MF = np.insert(MF, 0, Lectin[:, mkey], axis=1)
		
		
	return orderBook, Book
		
@register
def _atexit():
	print 'all DONE at: {}'.format(ctime())
		
if __name__ == '__main__':
	os.chdir("/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_CancerLectin/test/")
	
	'''
	# process CancerLectin and nonCancerLectin
	with open(r"./Dataset/cancerlectin178.txt", 'r') as f:
		seqs = [line.strip() for line in f if line[0] != '>']

	for seq in seqs:
		t, v = Tripeptide2Freq(seq)
		Lectin.append(v)
		
	with open(r"./Dataset/non-cancerlectin226.txt", 'r') as f:
		seqs = [line.strip() for line in f if line[0] != '>']
		
	for seq in seqs:
		t, v = Tripeptide2Freq(seq)
		Lectin.append(v)
	
	Lectin = np.array(Lectin, dtype=float)
	#svmfile(Lectin, './Dataset/data.txt', m, n)

	with open(r'./Dataset/Lectin', 'w') as f:
		pickle.dump(Lectin, f)
	'''	
	
	with open(r'./Dataset/Lectin', 'r') as f:
		Lectin = pickle.load(f)
		
	Lectin = np.array(Lectin, dtype=float)
	with open(r'result', 'r') as f:
		orderBook, Book = pickle.load(f)
	
	MF = np.empty(shape=(404, 1), dtype=float)
	for idx in orderBook:
		MF = np.insert(MF, 0, Lectin[:, idx], axis=1)
		
	orderBook, Book = optimalFunc(Lectin, Book, orderBook, MF, 4)
			
	
