#!/etc/bin/env python
# -*- coding: utf-8 -*-
#
# Initialization Dipeptide and Tripeptide component

import pickle

AminoAcid = ['G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','D','E','K','R','H']
Dipeptide = []
Tripeptide = []
for i in range(20):
	for j in range(20):
		Dipeptide.append(AminoAcid[i]+AminoAcid[j])
for i in range(20):
	for j in range(20):
		for k in range(20):
			Tripeptide.append(AminoAcid[i]+AminoAcid[j]+AminoAcid[k])

with open('../Dataset/Dipeptide', 'w') as f:
	pickle.dump(Dipeptide, f)

with open('../Dataset/Tripeptide', 'w') as f:
	pickle.dump(Tripeptide, f)
