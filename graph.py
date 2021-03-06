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

            rj = jj.findall(line)
            ri = ii.findall(line)

            mij = ij.findall(line)
            mji = ji.findall(line)

            #need to get p value from previous connection
            if(ri and rj):
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
                e = el.get_edge(sink, source, con_type)
                if(not e):
                    print(sink)
                    print(source)
                    print(con_type)
                    print("Earlier gap junction not found, exiting...")
                    exit(1)

            elif((mij and not mji) or (mji and not mij)):
                print("Repeat detected, but not for both cases?")

            else:
                #no special cases detected
                d = 1
                p_matches = p.findall(line)
                if(p_matches):
                    d = dec.findall(p_matches[0])[0]
                    d = float(d)

                el.new_edge(source, sink, connection_type = con_type, p = d)


    dotfile = open(dot_filename, 'w+')
    dotfile.write('digraph network {\n')


    cortex_set = (set(), set())
    thalamus_set = (set(), set())
    other_set = set()

    for node in node_set:
        found = False
        for item in cortex_options:
            if(item in node):
                if(node[-1] == 'L'):
                    cortex_set[0].add(node)
                elif(node[-1] == 'R'):
                    cortex_set[1].add(node)
                else:
                    print("Cortical neuron not lateralized? (%s) Exiting..." % node)
                    exit(1)
                found = True
                break
        if not found:
          for item in thalamus_options:
            if(item in node):
                if(node[-1] == 'L'):
                    thalamus_set[0].add(node)
                elif(node[-1] == 'R'):
                    thalamus_set[1].add(node)
                else:
                    print("Thalamic neuron not lateralized? (%s) Exiting..." % node)
                    exit(1)
                found = True
                break

        if not found:
            other_set.add(node)

    dotfile.write("subgraph cluster_cortex_left\n{\n    label = \"Cortex\"\n")
    dotfile.write("    style = filled\n")
    dotfile.write("    color = lightpink1\n")
    for node in cortex_set[0]:
        write_node(dotfile,node)
    dotfile.write("}\n")

    dotfile.write("subgraph cluster_cortex_right\n{\n    label = \"Cortex\"\n")
    dotfile.write("    style = filled\n")
    dotfile.write("    color = lightpink1\n")
    for node in cortex_set[1]:
        write_node(dotfile,node)
    dotfile.write("}\n")

    dotfile.write("subgraph cluster_thalamus_left\n{\n    label = \"Thalamus\"\n")
    dotfile.write("    style = filled\n")
    dotfile.write("    color = lightblue\n")
    for node in thalamus_set[0]:
        write_node(dotfile,node)
    dotfile.write("}\n")

    dotfile.write("subgraph cluster_thalamus_right\n{\n    label = \"Thalamus\"\n")
    dotfile.write("    style = filled\n")
    dotfile.write("    color = lightblue\n")
    for node in thalamus_set[1]:
        write_node(dotfile,node)
    dotfile.write("}\n")

    print("\nNodes not part of any predefined group:")
    for node in other_set:
        write_node(dotfile,node)
        print(node)



    el.write_edges(dotfile, cortex_set, thalamus_set)
    dotfile.write('}')
    f.close()
    dotfile.close()
    print("\nOutput written to: %s\n" % dot_filename)

def write_node(fp, node):
    fp.write("    %s [height=2, width=2, fontsize = 30, style = \"filled,solid\", fillcolor = white];\n" % node)



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
            print("Exiting")
            exit(1)

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

    def write_edges(self, fp, cortex_set, thalamus_set):
        for edge in self.edges:
            fp.write("    "  + edge.to_string(cortex_set, thalamus_set) + '\n')

class Edge:
    def __init__(self, source, sink, connection_type, p, fontsize):
        self.source = source
        self.sink = sink
        self.connection_type = connection_type
        self.p = p
        if(self.connection_type == 'gj'):
            self.bi = True
        else:
            self.bi = False
        self.fontsize = fontsize

    def to_string(self, cortex_set, thalamus_set):
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

            cortex_source = False
            cortex_sink = False
            for i in cortex_set:
                for item in i:
                    if self.source in item:
                        cortex_source = True
                        break
                    if self.sink in item:
                        cortex_sink = True
                        break

            if(cortex_source or cortex_sink):

                thalamus_source = False
                thalamus_sink = False

                for i in thalamus_set:
                    source_flag = False
                    sink_flag = False
                    for item in i:
                        if self.source in item:
                            thalamus_source = True
                            break
                        if self.sink in item:
                            thalamus_sink = True
                            break

            if(cortex_source and thalamus_sink):
                s = s + "color=red,"

            if(cortex_sink and thalamus_source):
                s = s + "color=blue,"

            s = s + "]"

        s = s + ";"
        return s

if __name__ == "__main__":
    main()
