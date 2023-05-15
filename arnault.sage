# https://core.ac.uk/download/pdf/81930829.pdf

def miller_rabin(n, base):
    assert n >= 3
    basis = list(primes(base))
    if n % 2 == 0:
        return False
    r, q = 0, n - 1
    while q % 2 == 0:
        q //= 2
        r += 1
    for base in basis:
        flag = True
        cur = pow(base, q, n)
        if cur == 1 or cur == n - 1:
            flag = False
        for _ in range(r-1):
            cur = pow(cur, 2, n)
            if cur == n - 1:
                flag = False
        if flag:
            return False
    return True

def getPseudoprime(base, bitLen, numTries):

	# Generate a pseudoprime with 3 factors that passes Miller-Rabin test up to base = 300

	baseList = list(primes(base)) 
	modsList = [4*prime for prime in baseList]

	# Note: this script only generates psudoprimes that are the product of 3 primes
	k1 = 1; k2 = 313; k3 = 353
	kValues = [k2,k3]
	kInverses = [k2-inverse_mod(k3,k2), k3-inverse_mod(k2,k3)]

	S = {}

	for i in range(len(baseList)):
		S[baseList[i]] = []	
		for j in range(3,modsList[i],2):
			if jacobi_symbol(baseList[i],j) == -1:
				S[baseList[i]].append(j)

	for i in range(len(baseList)):
		curMod = modsList[i]
		k1_inv = inverse_mod(k1,curMod); k2_inv = inverse_mod(k2,curMod); k3_inv = inverse_mod(k3,curMod)
		tmp1 = []; tmp2 = []; tmp3 = []
		for number in S[baseList[i]]:
			tmp1.append((k1_inv*(number+k1-1))%curMod)
			tmp2.append((k2_inv*(number+k2-1))%curMod)
			tmp3.append((k3_inv*(number+k3-1))%curMod)
		S[baseList[i]] = list(set(tmp1).intersection(set(tmp2).intersection(set(tmp3))))

	kValues = [k2,k3]
	kInverses = [k2-inverse_mod(k3,k2), k3-inverse_mod(k2,k3)]

	mod = lcm(modsList + kValues)
	while True:
		ZaList = [choice(S[prime]) for prime in baseList]
		curTry = crt(ZaList+kInverses, baseList+kValues) 

		for i in range(numTries):
			if i%(numTries//3) == 0:
				print ("Still searching...")
			p1 = getrandbits(bitLen) * mod + curTry
			p2 = k2*(p1 - 1) + 1
			p3 = k3*(p1 - 1) + 1
			n = p1*p2*p3
			if miller_rabin(n,base):
				return n,p1,p2,p3


def main():
	n,p1,p2,p3 = getPseudoprime(64,160,30000)
	print ("Found strong pseudoprime: ", n)
	print ("Factors: ", "\n--> p1: ", p1, "\n--> p2: ", p2, "\n--> p3: ", p3)
if __name__ == "__main__":
	main()
