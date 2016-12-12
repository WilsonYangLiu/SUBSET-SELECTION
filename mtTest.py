from atexit import register
from threading import Thread, Lock, current_thread
from time import ctime
from subprocess import call

class cleanOutputSet(set):
	def __str__(self):
		return ', '.join(x for x in self)

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

lock = Lock()
remaining = cleanOutputSet()

def svmfunc(fold):
	myname = current_thread().name
	with lock:
		remaining.add(myname)
		print '[{}] Started {}'.format(ctime(), myname)

	call(('grid.py -v {} ./heart.scale > ./tmp/{}'.format(fold, myname)), shell=True)
	freq = gridfile('./tmp/{}'.format(myname))

	with lock:
		remaining.remove(myname)
		print '[{}] Comleted {}'.format(ctime(), myname)
		print '    {}'.format(freq)
		print '    (remaining: {})'.format(remaining or None)

	return freq

def _main():
	lockFreq = Lock()
	Freq = {}
	for i in [4, 10]:
		freq = Thread(target=svmfunc, args=(i, )).start()
		with lockFreq:
			Freq[i] = freq


@register
def _atexit():
	print 'all DONE at: {}'.format(ctime())

if __name__ == '__main__':
	_main()
