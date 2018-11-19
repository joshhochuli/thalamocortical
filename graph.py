import re

def main():

    #enumerate membership for nodes, can't get from file
    cortex_options = ['FS', 'PY']
    thalamus_options = ['HTC', 'RE', 'RTC', 'IN']

    filename = 'thalcort_unilateral.py'
    f = open(filename)

    # find matches like "p = 0.8" 
    p = re.compile("p ?= ?[\',\"][0-9]?(\.[0-9]+)?[\',\"]")

    # find matches of 0.8, 2, etc
    dec = re.compile('[0-9]?(\.[0-9]+)')

    #find replication of earlier connection
    ii = re.compile('i ?= ?[^,]*i\[:\]')
    jj = re.compile('j ?= ?[^,]*j\[:\]')

    #find mirror of earlier connection
    ij = re.compile('i ?= ?[^,]*j\[:\]')
    ji = re.compile('j ?= ?[^,]*i\[:\]')

    dot_filename = 'test.dot'


    node_set = set()
    el = Edge_List()
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

            #need to get p value from previous connection
            if(ri and rj):
                print("Repeat earlier")
                el.repeat_edge(source,sink, con_type)

            #should never occur
            elif((ri and not rj) or (rj and not ri)):
                print("Repeat detected, but not for both cases?")
                print("ri:", ri)
                print("rj:", rj)
                print("exiting")
                exit(1)

            #should only occur if a previous gap juntion is defined, make sure
            #to check that this is the case
            elif(mij and mji):
                print("Mirror earlier")
                e = el.get_edge(sink, source, con_type)
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
                el.new_edge(source, sink, connection_type = con_type, p = d)


            curr_output_line = curr_output_line + ';\n'
            print(curr_output_line)
            print("--------")

    dotfile = open(dot_filename, 'w+')
    dotfile.write('digraph network {\n')


    cortex_set = set()
    thalamus_set = set()
    other_set = set()

    for node in node_set:
        if(node in cortex_options):
            cortex_set.add(node)
        elif(node in thalamus_options):
            thalamus_set.add(node)
        else:
            other_set.add(node)

    dotfile.write("subgraph cluster_cortex\n{\n")
    for node in cortex_set:
        write_node(dotfile,node)
    dotfile.write("}\n")

    dotfile.write("subgraph cluster_thalamus\n{\n")
    for node in thalamus_set:
        write_node(dotfile,node)
    dotfile.write("}\n")

    print("Nodes not part of any predefined group:")
    for node in other_set:
        write_node(dotfile,node)
        print(node)


    #dotfile.write("}\n")



    el.write_edges(dotfile)
    dotfile.write('}')
    f.close()
    dotfile.close()

def write_node(fp, node):
    fp.write("    %s [height=2, width=2, fontsize = 30];\n" % node)



class Edge_List:
    def __init__(self):
        self.edges = []

    def new_edge(self, source, sink, connection_type = None, p = 1, fontsize = 8):
        self.edges.append(Edge(source, sink, connection_type, p, fontsize))

    def mirror_edge(source, sink, connection_type):
        e = self.get_edge(source, sink, connection_type)

        if not e:
            print("edge not found")
            print(source)
            print(sink)
            print(connection_type)
            print(edges)

        self.new_edge(sink, source, connection_type, p = e.p, fontsize =
                e.fontsize)

    #mainly for finding and replicating earlier p value
    def repeat_edge(self, source, sink, connection_type):
        e = self.get_edge(source, sink)
     
        if not e:
            print("repeat edge not found")
            print(source)
            print(sink)
            print(connection_type)
            print(edges)
            print("exiting")
            exit(1)

        self.new_edge(source, sink, connection_type, p = e.p, fontsize =
                e.fontsize)

    def get_edge(self, source, sink, connection_type = None):
        if(connection_type):
            for edge in self.edges:
                if(edge.source == source and edge.sink == sink and edge.connection_type == connection_type):
                    return edge

            return None
        else:
            for edge in self.edges:
                if(edge.source == source and edge.sink == sink):
                    return edge

            return None

    def write_edges(self, fp):
        for edge in self.edges:
            fp.write("    "  + str(edge) + '\n')

class Edge:
    def __init__(self, source, sink, connection_type, p, fontsize):
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
if __name__ == "__main__":
    main()
