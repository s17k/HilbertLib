import random
random.seed()
n = input()
r = 0.1;
print r;
print n


def sting() :
	if random.randint(0,1) == 0 :
		return 1.0
	else :
 		return -1.0
def generuj() :
	x,y,z = random.random()*r*sting(), random.random()*r*sting(), random.random()*r*sting() 
	if x**2 + y**2 + z**2 <= r**2 :
		return x,y,z
	else :
		return generuj()

def rand_pos() :
	x,y,z = random.random(), random.random(), random.random()
	x *= sting()
	y *= sting()
	z *= sting()
	if x**2 + y**2 + z**2 <= (1-r)**2 :
		return x,y,z
	else :
		return rand_pos()

def zrob_kule(m) :
	a,b,c = rand_pos()
	while m > 0 :
		m-=1
		x,y,z = generuj()
		print x+a,y+b,z+c

for i in range(n) :
	print random.random()*sting(),random.random()*sting()

