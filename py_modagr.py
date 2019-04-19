#!/usr/bin/env python
#Manipulate a xmgrace .agr file, very simple and with basic operations only

from argparse import ArgumentParser 
import re

#This is the good color choice in order, the bright colors are not preferred as they are hard to read on the screen
list_xmgrace_color = [1, 2, 15, 4, 11, 6, 12, 13, 14, 10, 8, 3, 5, 7, 9]

#Graph properties
ar_graph_property = [
     "world",
     "stack world",
     "znorm",
     "view",
     "title",
     "title font",
     "title size",
     "title color",
     "subtitle",
     "subtitle font",
     "subtitle size",
     "subtitle color",
     "xaxes scale",
     "yaxes scale",
     "xaxes invert",
     "yaxes invert",
     "altxaxis",
     "altyaxis",
     "legend",
     "legend loctype",
     "legend",
     "legend box color",
     "legend box pattern",
     "legend box linewidth",
     "legend box linestyle",
     "legend box fill color",
     "legend box fill pattern",
     "legend font",
     "legend char size",
     "legend color",
     "legend length",
     "legend vgap",
     "legend hgap",
     "legend invert",
     "frame type",
     "frame linestyle",
     "frame linewidth",
     "frame color",
     "frame pattern",
     "frame background color",
     "frame background pattern",
     ]

ar_axis_property = [
     "axis",
     "axis type zero",
     "axis offset",
     "axis bar",
     "axis bar color",
     "axis bar linestyle",
     "axis bar linewidth",
     "axis label",
     "axis label layout",
     "axis label place",
     "axis label place",
     "axis label char size",
     "axis label font",
     "axis label color",
     "axis label place",
     "axis tick",
     "axis tick major",
     "axis tick minor ticks",
     "axis tick default",
     "axis tick place rounded",
     "axis tick",
     "axis tick major size",
     "axis tick major color",
     "axis tick major linewidth",
     "axis tick major linestyle",
     "axis tick major grid",
     "axis tick minor color",
     "axis tick minor linewidth",
     "axis tick minor linestyle",
     "axis tick minor grid",
     "axis tick minor size",
     "axis ticklabel",
     "axis ticklabel format",
     "axis ticklabel prec",
     "axis ticklabel formula",
     "axis ticklabel append",
     "axis ticklabel prepend",
     "axis ticklabel angle",
     "axis ticklabel skip",
     "axis ticklabel stagger",
     "axis ticklabel place",
     "axis ticklabel offset",
     "axis ticklabel start type",
     "axis ticklabel start",
     "axis ticklabel stop type",
     "axis ticklabel stop",
     "axis ticklabel char size",
     "axis ticklabel font",
     "axis ticklabel color",
     "axis tick place",
     "axis tick spec type",
     "axis tick spec",
]

ar_set_property = [
        "hidden",
        "type",
        "symbol",
        "symbol size",
        "symbol color",
        "symbol pattern",
        "symbol fill color",
        "symbol fill pattern",
        "symbol linewidth",
        "symbol linestyle",
        "symbol char",
        "symbol char font",
        "symbol skip",
        "line type",
        "line linestyle",
        "line linewidth",
        "line color",
        "line pattern",
        "baseline type",
        "baseline",
        "dropline",
        "fill type",
        "fill rule",
        "fill color",
        "fill pattern",
        "avalue",
        "avalue type",
        "avalue char size",
        "avalue font",
        "avalue color",
        "avalue rot",
        "avalue format",
        "avalue prec",
        "avalue prepend",
        "avalue append",
        "avalue offset",
        "errorbar",
        "errorbar place",
        "errorbar color",
        "errorbar pattern",
        "errorbar size",
        "errorbar linewidth",
        "errorbar linestyle",
        "errorbar riser linewidth",
        "errorbar riser linestyle",
        "errorbar riser clip",
        "errorbar riser clip length",
        "comment",
        "legend"
        ]

#Join axis property to graph property
for axis in ("x", "y"):
    ar_graph_property += [axis + x for x in ar_axis_property]

#Check properties are both graph and set, and forbids them in set_property_auto
ar_both_property = [x for x in ar_graph_property if x in ar_set_property]

class agrLine(object):
    '''
    A object represent one line in agr file
    '''
    re_set = re.compile(r"s(\d+)")

    def __init__(self, line=None):
        self.s = None
        self.g = None
        self.name = None
        self.value = None
        self.command = None #A command without arguments
        if (line is not None):
            self.parse(line)
        return

    @property
    def is_graph(self):
        '''
        Is a property of a graph (not a set)
        '''
        return self.s is None

    @property
    def is_set(self):
        '''
        Is a property of a set (not a graph)
        '''
        return self.s is not None

    def parse(self, line):
        '''
        Parse a line in an AGR file, split into name and value
        '''
        line = line.strip()
        #print("parse -- ", line)
#Not start with @, not a command
        if (line[0] != "@"):
            return
        l = line[1:].replace(",", " , ")
        ar = []

        ix_string_start = -1
        ix_term_start = -1

#We do not use split as we may have a string
        l = l + " " #Add a space to end the last term
        for i in range(len(l)):
            if (l[i] == '"'):#String
                if (ix_string_start >= 0):
                    ar.append(l[ix_string_start:i+1])
                    ix_string_start = -1
                    ix_term_start = -1 #Should no term also
                else:
                    ix_string_start = i
            elif (l[i] == " "):
                if (ix_string_start >= 0):
                    continue
                elif (ix_term_start >= 0):
                    ar.append(l[ix_term_start:i])
                    ix_term_start = -1
            elif (ix_term_start < 0):
                ix_term_start = i
            else:
                continue
        #print(ar)

#Find comma splitted value and join them
        n_firstcomma = None
        for i, x in enumerate(ar):
            if ("," in x):
                n_firstcomma = i
                break

#Find set index and remove "s*" from name
        m = agrLine.re_set.match(ar[0])
        if (m is not None):
            self.s = int(m.group(1))
            ar = ar[1:]

        if (n_firstcomma is not None):
            self.name = " ".join(ar[:n_firstcomma-1])
            self.value = [x for x in ar[n_firstcomma-1:] if x != ","]
        elif (len(ar) > 1):
            self.name = " ".join(ar[:-1])
            self.value = ar[-1]
        else:#Only one term, it is a command
            self.command = ar[-1]

        return

    def __str__(self):
        '''
        Convert to line
        '''
        if (self.command is not None):
            return "@    %s" % self.command
        else:
            if (isinstance(self.value, list)):
                v = ", ".join(self.value)
            else:
                v = self.value

            if (self.is_set):
                k = "s%i %s" % (self.s, self.name)
            else:
                k = self.name
            
            return "@    %s %s" % (k, v)

def split_multiline_to_single(lines):
    '''
    Split an array with multi-line strings to an array with single line string
    This is useful when we want to process line by line
    '''
#lines may contain multi-line string, split them
    lines2 = []
    for line in lines:
        if (line[-1] != "\n"):
            line = line + "\n"

        i0 = 0
        i = line.find("\n", i0)
        while (i < len(line) - 1 and i != -1):
            i += 1
            lines2.append(line[i0:i])
            i0 = i
            i = line.find("\n", i0)

        lines2.append(line[i0:])

    return lines2


def set_axis_limit(lines, xmin=None, ymin=None, xmax=None, ymax=None, graph=None):
    '''
    Set the axis limitations

    :param lines: may be modified in-place; use return value as the result
    :return: the modified lines
    '''
    if (xmin is None and ymin is None and xmax is None and ymax is None):
        return lines
    world_new = [xmin, ymin, xmax, ymax]
    for i in range(len(lines)):
        line = lines[i]
        if ("with g" in line):
            n_graph = float(line[line.rfind("g")+1:])
        if ("world " in line and (graph is None or n_graph == graph)):
            world = [float(x) for x in line[line.find("world")+5:].split(",") if x.strip() != ""]
            world = [vnew if vnew is not None else v for (v, vnew) in zip(world, world_new)]
            lines[i] = "@    world " + ", ".join([str(x) for x in world]) + "\n"

    return lines

def set_property_auto(lines, dic_property, ix_graph=None, ix_set=None):
    '''
    Set properties, automatically treat different key in SET or GRAPH
    Also support xmin, xmax ymin, ymax as graph properties
    To set s* properties, at least one s* command must present in the original lines
    '''
    dic_property_set = {}
    dic_property_graph = {}
    dic_world = {}

    ar = ["xmin", "ymin", "xmax", "ymax"]

    for key, val in dic_property.items():
        if (key in ar_both_property):
            raise ValueError("Property %s is ambiguous for graph and set, use set_property instead")
        if (key in ar_set_property):
            dic_property_set[key] = val
        elif (key in ar_graph_property):
            dic_property_graph[key] = val
        elif (key in ar):
            dic_world[key] = val
        else:
            print("Skip unknown property %s" % key)
#           raise ValueError("Unknown property \"%s\"" % key)

    lines = set_axis_limit(lines, graph=ix_graph, **dic_world)
    lines = set_property(lines, dic_property_graph, ix_graph, dic_property_set, ix_set)

    return lines



def set_property(lines, dic_property_graph={}, ix_graph=None, dic_property_set={}, ix_set=None):
    '''
    Set properties of a given graph or a given set
    To set s* properties, at least one s* command must present in the original lines

    :param lines: may be modified in-place; use return value as the result
    :param dic_property_graph: a dictionary of property names and values. If the value is None, then the original one is kept
    :return: the modified lines
    '''
    if (len(dic_property_graph) == 0 and len(dic_property_set) == 0):
        return lines

    def mark_dic_found(dic1):
        '''
        Convert a dictionary into all False dictionary
        None is ignored
        '''
        dic2 = {}
        for key, val in dic1.items():
            if (val is not None):
                dic2[key] = False
        return dic2

    def insert_lines_from_dic(dic_property, dic_found, lines, n_startline, func_format):
        '''

        :param i: Insert after line i
        :return: new lines and how many lines are inserted 
        '''
        lines2 = []
        for key, val in dic_found.items():
            if (not val):
                lines2.append(func_format(key, dic_property[key]))
#                lines2.append("@    %s %s\n" % (key, dic_property[key]))
        if (len(lines2) != 0):
            lines = lines[:n_startline] + lines2 + lines[n_startline:]

        return lines, len(lines2)


    i = 0
    dic_found_graph = {}
    dic_found_set = {}
    n_graph = -1
    n_set = -1
    re_g_on = re.compile("@g(\\d+) on\n")
    while (i <= len(lines)):
        if (i == len(lines)):
            line = ""
        else:
            line = lines[i]
#End of the last graph (new graph, or EOF reached)
#Note, sometimes new graph start with "with g*" and sometimes "g* on" then "with g*"
        m1 = re_g_on.match(line)
#If all lines are read or meet @target, we have done
        if ("with g" in line or m1 is not None or i == len(lines) or "@target" in line):
            n_graph_read = -2
            b_end = False
            if ("with g" in line):
                n_graph_read = int(line[line.rfind("g")+1:])
            elif (m1 is not None):
                n_graph_read = int(m1.group(1))
            else: #Meet EOF or @target, all options finish
                b_end = True
            if (n_graph_read == n_graph): #Already started, just skip
                i += 1
                continue
#Check if all items are set, if not, set them for the previous set
#For both graph and set
            if (n_set >= 0 and (ix_set is None or ix_set == n_set)):
                lines, i2 = insert_lines_from_dic(dic_property_set, dic_found_set, lines, i,
                    lambda key, val:"@    s%i %s %s\n" % (n_set, key, val))
                i += i2
            if (n_graph >= 0):
                lines, i2 = insert_lines_from_dic(dic_property_graph, dic_found_graph, lines, n_startline,
                        lambda key, val:"@    %s %s\n" % (key, val))
                i += i2
#Set up for a new graph
            if (b_end):
                break 

#Process new graph
            n_graph = n_graph_read
            n_startline = i + 1

#Use a dictionary to record whether an item is found or not
            dic_found_graph = mark_dic_found(dic_property_graph)
            dic_found_set = mark_dic_found(dic_property_set)
        else:
            if (ix_graph is None or ix_graph == n_graph):
                la = agrLine(line)
#Set property end
                if (la.s is not None and la.s != n_set ):
                    if (n_set >= 0 and (ix_set is None or ix_set == n_set)):
                        lines, i2 = insert_lines_from_dic(dic_property_set, dic_found_set, lines, i,
                            lambda key, val:"@    s%i %s %s\n" % (n_set, key, val))
                        i += i2
                    n_set = la.s
#Graph property
                if (la.s is None and dic_property_graph.has_key(la.name) 
                        and (n_graph == ix_graph or ix_graph is None)):
                    la.value = dic_property_graph[la.name]
                    lines[i] = str(la) + "\n"
                    dic_found_graph[la.name] = True
#Set property
                elif (la.s is not None and dic_property_set.has_key(la.name) and 
                        (n_set == ix_set or ix_set is None)):
                    la.value = dic_property_set[la.name]
                    lines[i] = str(la) + "\n"
                    dic_found_set[la.name] = True

        i += 1

    return lines

def sprintf_agr(list_data, list_legend=None, b_col=None):
    '''
    Convert a list of xy data to xmgrace lines

    :param list_legend: the legend to plot
    :param b_col: data in list_data is col-wise or row-wise
    '''
    lines = []
    lines.append("@g0 on\n")
    lines.append("@with g0\n")
    if (list_legend is not None):
        for i, legend in enumerate(list_legend):
            lines.append("@    s%i legend \"%s\"\n" % (i, legend))
    for i, data in enumerate(list_data):
        lines.append("@target G0.S%i\n" % i)
        lines.append("@type xy\n")
        if (b_col):
            for x,y in zip(*data):
                lines.append("%f %f\n" % (x,y))
        else:
            for x,y in data:
                lines.append("%f %f\n" % (x,y))
        lines.append("&\n")
    return lines


def plot_agr_with_legend(filename, list_data, list_legend, b_col = True):
    '''
    Plot data with legend names to given file

    :param filename: the output filename
    :param list_data: a list of xy data, in col or row mode
    :param list_legend: legend name for each data
    :param b_col: where data in list_data is col wise or row wise
    '''
    with open(filename, 'w') as f:
        f.write("@g0 on\n")
        f.write("@with g0\n")
        for i, legend in enumerate(list_legend):
            f.write("@    s%i legend \"%s\"\n" % (i, legend))
        for i, data in enumerate(list_data):
            f.write("@target G0.S%i\n@type xy\n" % i)
            if (b_col):
                for x,y in zip(*data):
                    f.write("%f %f\n" % (x,y))
            else:
                for x,y in data:
                    f.write("%f %f\n" % (x,y))
            f.write("&\n")

    return


def Main():
    description = '''This script does minimum modification to a xmgrace .agr file by editing, adding or deleting some lines.
Note to make it work, there must be at least one s* command in agr file to make modification of s* to work, which is always guaranteed if this file is produced by xmgrace.
Pay attention to that the original file is overwrite!'''
    parser = ArgumentParser(description=description)

    parser.add_argument("filenames", type=str, nargs="*", help="Filenames that will be modified")
    #, with -i, this can modifies multiple filenames, without -i, there must be two filenames, read from the first one and save to the second one")
    #parser.add_argument("-i", "--inplace", action="store_true", default=False, dest="b_inplace", help="If set, modify filename inplace instead rename it ")
    parser.add_argument("-g", type=int, help="The index of the graph (start from 0)")
    parser.add_argument("-s", type=int, help="The index of the graph (start from 0)")
    parser.add_argument("--xmin", type=float, help="Minimum x axis")
    parser.add_argument("--xmax", type=float, help="Maximum x axis")
    parser.add_argument("--ymin", type=float, help="Minimum y axis")
    parser.add_argument("--ymax", type=float, help="Maximum y axis")
#Below are some special keys
    parser.add_argument("--all-linewidth", type=float, help="set frame-linewidth, line-linewidth and symbol-linewidth")
    parser.add_argument("--frame-linewidth", type=float, help="Linewidth of the frame")
    parser.add_argument("--line-linewidth", type=float, help="Linewidth of the line of a set")
    parser.add_argument("--symbol-linewidth", type=float, help="Linewidth of the symbol")
    parser.add_argument("--symbol-size", type=float, help="Symbol size of the set")

#Some general keys
    parser.add_argument("-p", "--property", type=str, help='A filename for the format (without space) or a string in the format of "key value", seperated with ";", for example --property "xaxis tick on;line color 1". The file format is the same, but each line contains a "key value".')

    options = parser.parse_args()
#Set all line width
    if (options.all_linewidth is not None):
        options.frame_linewidth = options.all_linewidth
        options.line_linewidth = options.all_linewidth
        options.symbol_linewidth = options.all_linewidth

    def conv_option(ar):
        '''
        Convert "_" in options to " "
        '''
        dic = {}
        for x in ar:
            v = getattr(options ,x)
            if (v is not None):
                dic[x.replace("_", " ")] = v
        return dic
                

    dic_property_graph = conv_option(["frame_linewidth"])
    dic_property_set = conv_option(["symbol_linewidth", "line_linewidth", "symbol_size"])

#Read --property
    if (options.property is not None):
        if (" " in options.property): #A string
            ar1 = [x for x in options.property.split(";") if x != ""]
        else:
            with open(options.property, "r") as f:
                ar1 = [x.strip() for x in f.readlines()]
                ar1 = [x for x in ar1 if x != ""]
        
        def split_value(st):
            '''
            Convert a "name value" with spaces in name to (name, val)
            '''
            st = st.strip()
            i = st.rfind(" ")
            return (st[:i], st[i+1:])

        for x in ar1:
            key, val = split_value(x)
            if (key in ar_set_property):
                dic_property_set[key] = val
            elif (key in ar_graph_property):
                dic_property_graph[key] = val
            else:
                raise ValueError("Unknown property \"%s\"" % key)


    for filename in options.filenames:
        #print("Start parse %s" % filename)
        with open(filename, 'r') as f:
            lines = f.readlines()

        lines = set_axis_limit(lines, options.xmin, options.ymin, options.xmax, options.ymax, options.g)

        lines = set_property(lines, dic_property_graph, options.g, dic_property_set, options.s)
        
#Write into file
        with open(filename ,'w') as f:
            f.writelines(lines)


if __name__ == "__main__":
    Main();
