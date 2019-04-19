#!/usr/bin/env python
from argparse import ArgumentParser
import re


parser = ArgumentParser(description="Combine several .agr files and distinguish datasets in different files with line color")
parser.add_argument('filename',type=str,nargs="*",help="AGR files")
parser.add_argument("-o",dest='outputfilename',type=str,required=False,default="comb.agr",help="Output filename")
args = parser.parse_args()

ix_add = 0
re1 = re.compile("([ \.]s)(\d+)")
re_color = re.compile("([ \.]s\d+ line color )(\d+)")
s_max = 0

list_color = [1, 2, 4, 6, 8, 10, 11, 12, 13, 14, 15]

def add_s(m):
    global s_max
    s = int(m.group(2))
    if (s> s_max):
        s_max = s
    return m.group(1) + str(s+ix_add)

def mod_color(m):
    global ix_file
    return m.group(1) + str(list_color[ix_file])

with open(args.outputfilename, 'w') as f:
    for ix_file, filename in enumerate(args.filename):
        s_max = 0
        with open(filename, 'r')  as f2:
            lines = f2.readlines()
#Replace all sx to sx+ix
        for i in xrange(len(lines)):
            lines[i] = re1.sub(add_s,lines[i])
            lines[i] = re_color.sub(mod_color,lines[i])
        
        ix_add += s_max + 1

        for line in lines:
            f.write(line)


print("Merge files (%s) into %s " % (" ".join(args.filename), args.outputfilename))



