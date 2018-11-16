import re

filename = 'thalcort_unilateral.py'
f = open(filename)

# find matches like "p = 0.8" 
p = re.compile("p ?= ?[\',\"][0-9]?(\.[0-9]+)?[\',\"]")

# find matches of 0.8, 2, etc
dec = re.compile('[0-9]?(\.[0-9]+)')

#special connection cases
#replicate earlier connection
ii = re.compile('i ?= ?[^,]*i\[:\]')
jj = re.compile('j ?= ?[^,]*j\[:\]')

#mirror earlier connection
ij = re.compile('i ?= ?[^,]*j\[:\]')
ji = re.compile('j ?= ?[^,]*i\[:\]')

for line in f:

    #only concerned with actual connections
    if('connect' in line):

        type = None
        print(line)
        s = line.split('.')
        t = s[0].split('_')
        sink = t[0]
        source = t[1]

        #type might not be specified, like with RND case
        if(len(t) == 3):
            type = t[2]

        print("Source: %s" % source)
        print("Sink: %s" % sink)

        if(type):
            print("Type: %s" % type)
        else:
            print("Type: random")

        p_matches = p.findall(line)
        if(p_matches):
            print("has a p")
            d = dec.findall(p_matches[0])[0]

            print("P: %.2f" % float(d))

        rj = jj.findall(line)
        ri = ii.findall(line)

        mij = ij.findall(line)
        mji = ji.findall(line)

        if(ri and rj):
            print("Repeat earlier")

        if((ri and not rj) or (rj and not ri)):
            print("Repeat detected, but not for both cases?")
            print("ri:", ri)
            print(rj)

        if(mij and mji):
            print("Mirror earlier")

        if((mij and not mji) or (mji and not mij)):
            print("Repeat detected, but not for both cases?")

        print("--------")

