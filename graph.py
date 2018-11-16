import re
edges = []

class Edge:
    def __init__(self, source, sink, connection_type = None, p = 1, fontsize = 8):
        global edges
        self.source = source
        self.sink = sink
        self.connection_type = connection_type
        self.p = p
        print(self.p)
        if(self.connection_type == 'gj'):
            self.bi = True
        else:
            self.bi = False
        self.fontsize = fontsize
        edges.append(self)

    def to_string(self):
        s = self.source + " -> " + self.sink
        if(self.bi):
            direction = "dir=\"both\","
        else:
            direction = ""
        if(self.p or self.fontsize or self.bi):
            s = s + " ["
            if(self.p):
                if(self.connection_type):
                    s = s + "label = \"%s, \\np = %.2f\",%s" % (self.connection_type, self.p, direction)
                else:
                    s = s + "label = \"\\np = %.2f\"," % (self.p)
            if(self.fontsize):
                s = s + "fontsize = %d," % self.fontsize
            if(self.source == "RND" or self.sink == "RND"):
                s = s + "style=dashed,"
            s = s + "]"

        s = s + ";"
        return s

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string()

def mirror_edge(source, sink, connection_type):
    global edges
    e = get_edge(source, sink, connection_type)

    if not e:
        print("edge not found")
        print(source)
        print(sink)
        print(connection_type)
        print(edges)

    Edge(sink, source, connection_type, p = e.p, fontsize =
            e.fontsize)

def repeat_edge(source, sink, connection_type):
    e = get_edge(source, sink)
 
    if not e:
        print("repeat edge not found")
        print(source)
        print(sink)
        print(connection_type)
        print(edges)

    Edge(source, sink, connection_type, p = e.p, fontsize =
            e.fontsize)

def get_edge(source, sink, connection_type = None):
    global edges
    if(connection_type):
        for edge in edges:
            if(edge.source == source and edge.sink == sink and edge.connection_type == connection_type):
                return edge

        return None
    else:
        for edge in edges:
            if(edge.source == source and edge.sink == sink):
                return edge

        return None


def write_edges(fp):
    global edges
    for edge in edges:
        fp.write(str(edge) + '\n')

def main():
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

    dot_filename = 'test.dot'

    dotfile = open(dot_filename, 'w+')
    dotfile.write('digraph network {\n')

    node_set = set()
    global edges

    #for clustering output into different sections
    #subgraph names must be start with "cluster" to be grouped
    #hard-coded, no way to determine clustering from file
    dotfile.write(""" subgraph cluster_cortex\n{\n  FS;\n  PY;\n  color=blue;\n}
     subgraph cluster_thalamus\n{\n  HTC;\n  RTC;\n  RE;\n  IN;\n}\n""")


    for line in f:

        #only concerned with actual connections
        if('connect' in line):

            curr_output_line = ''
            con_type = None
            print(line)
            s = line.split('.')
            t = s[0].split('_')
            sink = t[0]
            source = t[1]


            if sink not in node_set:
                #create node
                node_set.add(sink)

            if source not in node_set:
                #create node
                node_set.add(source)


            #type might not be specified, like with RND case
            if(len(t) == 3):
                con_type = t[2]

            print("Source: %s" % source)
            print("Sink: %s" % sink)

            if(con_type):
                print("Type: %s" % con_type)
            else:
                print("Type: random")
            rj = jj.findall(line)
            ri = ii.findall(line)

            mij = ij.findall(line)
            mji = ji.findall(line)

            if(ri and rj):
                print("Repeat earlier")
                repeat_edge(source,sink, con_type)

            elif((ri and not rj) or (rj and not ri)):
                print("Repeat detected, but not for both cases?")
                print("ri:", ri)
                print(rj)

            elif(mij and mji):
                print("Mirror earlier")
                e = get_edge(sink, source, con_type)
                if(not e):
                    print("Earlier gap junction not found, exiting...")
                    exit(1)

            elif((mij and not mji) or (mji and not mij)):
                print("Repeat detected, but not for both cases?")

            else:
                #no special cases detected
                #e = gv.edge(source, sink, source + sink)
                curr_output_line = curr_output_line + '%s -> %s' % (source, sink)

                d = 1
                p_matches = p.findall(line)
                if(p_matches):
                    print("has a p")
                    d = dec.findall(p_matches[0])[0]
                    d = float(d)

                    print("P: %.2f" % float(d))
                    curr_output_line = curr_output_line + ' [label = \"p = %s\"]' % d


                print("e source: %s" % source)
                print("ct: %s" % con_type)
                Edge(source, sink, con_type, p = d)


            curr_output_line = curr_output_line + ';\n'
            print(curr_output_line)
            print("--------")

    write_edges(dotfile)
    dotfile.write('}')
    f.close()
    dotfile.close()

    print(edges)

if __name__ == "__main__":
    main()
