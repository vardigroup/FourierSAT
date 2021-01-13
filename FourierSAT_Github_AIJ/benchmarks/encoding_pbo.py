import itertools
import sys

def cnf2pb(filepath,topath):
    f = open(filepath,'r+')
    g = open(topath,'w+')
    flines = f.readlines()
    for line in flines:
        split = line.split()
        if len(split) == 0: continue
        if split[0] == 'p':
            n = int(split[2])
            m = int(split[3])
            g.write("* #variable= "+repr(n)+" #constraint= " +repr(m)+' ;\n')
        elif split[0] == 'g':
            for i in range(n):
                g.write('+1 x'+repr(i+1)+' ')
            g.write('<= '+repr(abs(int(split[1])))+' ;\n')
        else:
            for i in range(len(split)-1):
                g.write('+1 x'+repr(int(split[i]))+' ')
            g.write('>= 1 ;\n')

    f.close()
    g.close()

cnf2pb(sys.argv[1],sys.argv[2])
