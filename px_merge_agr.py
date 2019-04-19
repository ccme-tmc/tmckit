#!/usr/bin/env python
from argparse import ArgumentParser
import re


parser = ArgumentParser(description="Combine several .agr files and distinguish datasets in different files with line color")
parser.add_argument('filename',type=str,nargs="*",help="AGR files")
parser.add_argument("-o",dest='outputfilename',type=str,required=False,default="comb.agr",help="Output filename")
parser.add_argument("--legend",dest='legendname',type=str,required=False,default=None,help="Legend names for each .agr, splitted with spaces")
parser.add_argument("--plottype",dest='plottype',type=str,required=False,default=None,help="Symbol or line names for each .agr, splitted with space, 1(line), 2(symbol), 3(both)")
args = parser.parse_args()

b_has_legends = args.legendname is not None
if (b_has_legends):
    legends = args.legendname.split()
else:
    legends = [None for x in args.filename]

b_has_plottype = args.plottype is not None
if (b_has_plottype):
    plottypes = [int(x) for x in args.plottype.split()]
    if (len(plottypes) != len(args.filename)):
        raise ValueError("Number of plottypes != filenames")
#   plottypes = [(plottypes[2*i], plottypes[2*i+1]) for i in xrange(len(plottypes)/2)]
    plottypes = [(x, i) for i,x in enumerate(plottypes)]
else:
    plottypes = [(1,i) for i,x in enumerate(args.filename)]
#Convert symbol as 1 and line as 2
#Group two as one

ix_add = 0
re1 = re.compile(r"([ \.]s)(\d+)")
re_all = re.compile(r"(@\s+s\d+) (.+) (\S+)")
re_color = re.compile(r"([ \.]s\d+ line color )(\d+)")
re_linetypecolor = re.compile(r"([ \.]s\d+ line color )(\d+)")
re_legend = re.compile(r"([ \.]s\d+ legend ).*")
s_max = 0

list_color = [1, 2, 4, 15, 6, 8, 10, 11, 12, 13, 14]

def change_line_symbol(plottype, color, line):
    '''
    Change lines to given plottype and color 
    '''
    match = re_all.match(line)
    if (match is None):
        return line
#       raise ValueError("Cannot parse given value %s" % line.strip())

    name = match.group(2)
    line2 = "%s %s " % (match.group(1), match.group(2))
    if (name == "symbol"):
        if (plottype == 2 or plottype == 3):
            line2 += "1"
        else:
            line2 += "0"
    elif (name == "line type"):
        if (plottype == 1 or plottype == 3):
            line2 += "1"
        else:
            line2 += "0"
    elif (name == "line color" or name == "symbol color"):
        line2 += str(list_color[color])
    elif (name == "symbol size"):
        line2 += str(0.4)

    return line2 + "\n"

def add_s(m):
    global s_max
    s = int(m.group(2))
    if (s> s_max):
        s_max = s
    return m.group(1) + str(s+ix_add)

def mod_color(m):
    global ix_file
    return m.group(1) + str(list_color[ix_file])

def mod_legend(m):
    global legend_now
    return m.group(1) + '"' + legend_now + '"'

with open(args.outputfilename, 'w') as f:
    for ix_file, (filename, legend_now) in enumerate(zip(args.filename, legends)):
        s_max = 0
        with open(filename, 'r')  as f2:
            lines = f2.readlines()
#Replace all sx to sx+ix
        for i in xrange(len(lines)):
            lines[i] = re1.sub(add_s,lines[i])
            if (b_has_legends):
                lines[i] = lines[i].replace("legend off", "legend on")
#Legend: replace only once
            if (legend_now is not None):
                m1 = re_legend.match(lines[i])
                if (m1 is not None):
                    lines[i] = re_legend.sub(re_legend,lines[i])
                    legend_now = None
            lines[i] = change_line_symbol(plottypes[ix_file][0], plottypes[ix_file][1], lines[i])
#           lines[i] = re_color.sub(mod_color,lines[i])
#Not replaced, add
        if (legend_now is not None):
            lines.append("@   s%i legend \"%s\"\n" % (ix_add, legend_now))
        
        ix_add += s_max + 1

        for line in lines:
            f.write(line)


print("Merge files (%s) into %s " % (" ".join(args.filename), args.outputfilename))



