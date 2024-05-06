import numpy as np
import copy

def countingNucleotides(fileName):
    f = open(fileName, "r")
    dna = f.readline()
    a = dna.count("A")
    c = dna.count("C")
    g = dna.count("G")
    t = dna.count("T")
    return "{} {} {} {}".format(a, c, g, t)

####################################################

def transcribeDNA2RNA(fileName):
    f = open(fileName, "r")
    dna = f.readline()
    res = open("rosalind_rna_res.txt", "w")
    res.write(dna.replace("T", "U"))
    return dna.replace("T", "U")

####################################################

def DNAcomplement(fileName):
    f = open(fileName, "r")
    dna = f.readline()

    res = ""
    basePairs = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    for base in dna[::-1]:
        res += basePairs.get(base)

    resF = open("test_res.txt", "w")
    resF.write(res)

    return res
    
####################################################

def fibWrapper(fileName):
    f = open(fileName, "r")
    n, k = f.readline().split(" ")

    res = fibonacciRabbits(int(n), int(k), 1, 0)

    resF = open("test_res.txt", "w")
    resF.write(str(res))

    return res

def fibonacciRabbits(n, k, baby, adult):
    if n == 1:
        return baby+adult
    else:
        return fibonacciRabbits(n-1, k, baby=adult*k, adult=adult+baby)

def mortalFibWrapper(fileName):
    f = open(fileName, "r")
    n, m = f.readline().split(" ")

    res = mortalFibRabbits(int(n), int(m), [[1,1]])

    resF = open("test_res.txt", "w")
    resF.write(str(res))

    return res

def mortalFibRabbits(n, m, lifespan):
    if n==1:
        return np.sum(lifespan, axis=0)[0]
    else:
        newLifespan = updateRabbits(lifespan, m)
        return mortalFibRabbits(n-1, m, newLifespan)

def updateRabbits(arr, m):
    res = []
    babyCount = 0
    for numRabbits, months in arr:
        if months < m: res.append([numRabbits, months+1])
        if months != 1: babyCount += numRabbits
    if babyCount > 0: res.append([babyCount, 1])
    return res

####################################################

def mendelInheritWrapper(fileName):
    f = open(fileName, "r")
    k, m, n = f.readline().split(" ")

    res= mendelianInheritance(int(k), int(m), int(n))

    resF = open("test_res.txt", "w")
    resF.write(str(res))

    return res

def mendelianInheritance(k, m, n):
    organisms = "k" * k + "m" * m + "n"*n
    orgLst = list(organisms)

    totalProb = 0
    
    prob = {
        "k": ['1','1'],
        "m": ['1','0'],
        "n": ['0','0']
    }

    for org in orgLst:
        modOrgLst = copy.deepcopy(orgLst)
        modOrgLst.remove(org)

        if org == "k":
            totalProb += 1 * len(modOrgLst)

        else: 
            for org2 in modOrgLst:
                if org2 == "k": totalProb += 1
                else:
                    alleles = prob.get(org)
                    alleles2 = prob.get(org2)
                    totalProb += checkAlleles(alleles, alleles2)
                
    return totalProb/(len(orgLst)*(len(orgLst)-1))

def checkAlleles(a1, a2):
    p = 0
    for a in a1:
        for b in a2:
            temp = a+b
            if "1" in temp: p += 1/4
    return p

####################################################

def expectedOffspringWrapper(fileName):
    f = open(fileName, "r")
    a, b, c, d, e, f = f.readline().split(" ")

    res = expectedOffspring(int(a), int(b), int(c), int(d), int(e), int(f))

    resF = open("test_res.txt", "w")
    resF.write(str(res))

    return res

def expectedOffspring(a, b, c, d, e, f):
    population = [a, b, c, d, e, f]

    k = len(population)
    res = 0
    for i in range(1,k+1):
        res += p(i, population[i-1]) * 2
    return res

def p(x, magX):
    probs = {
        1: 1.,
        2: 1.,
        3: 1.,
        4: 3/4,
        5: 1/2,
        6: 0
    }

    return probs.get(x) * magX

####################################################

