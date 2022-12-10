import copy
import time


class HyperEdge:

    def __init__(self, id, vertices, label, weight):
        self.id = id
        self.vertices = vertices
        self.label = label
        self.weight = weight

    def get_id(self):
        return self.id

    def get_vertices(self):
        return self.vertices

    def get_label(self):
        return self.label

    def equals(self, hyper_edge):
        if self.label != hyper_edge.label:
            return False
        if len(self.vertices) != len(hyper_edge.vertices):
            return False
        for v in self.vertices:
            if v not in hyper_edge.vertices:
                return False
        return True

    def __str__(self):
        result = dict()
        result["id"] = str(self.id)
        result["vertices"] = str(self.vertices)
        result["label"] = str(self.label)
        return str(result)


class HyperGraph:

    def __init__(self):
        self.vertex_label = dict()
        self.vertex_label_bucket = dict()
        self.edge_label_bucket = dict()
        self.hyper_edges = []
        self.vertex_count = 0
        self.edge_count = 0
        self.weight = 0

    def find_isomorphisms(self, H):
        v_isos, rev_v_isos, e_isos, rev_e_isos = [], [], [], []
        n_v_isos, n_rev_v_isos = self.find_isomorphisms_v(0, list(self.vertex_label.keys()), dict(), dict(), H)

        for i in range(len(n_v_isos)):
            v_iso = n_v_isos[i]
            rev_v_iso = n_rev_v_isos[i]
            n_e_isos, n_rev_e_isos = self.find_isomorphisms_e(0, self.hyper_edges, dict(), dict(), v_iso, rev_v_iso, H)
            for j in range(len(n_e_isos)):
                e_iso = n_e_isos[j]
                rev_e_iso = n_rev_e_isos[j]
                v_isos.append(copy.deepcopy(v_iso))
                rev_v_isos.append(copy.deepcopy(rev_v_iso))
                e_isos.append(e_iso)
                rev_e_isos.append(rev_e_iso)
        return v_isos, rev_v_isos, e_isos, rev_e_isos

    def find_isomorphisms_e(self, i, hyperedges, e_iso, rev_e_iso, v_iso, rev_v_iso, H):
        if i == len(hyperedges):
            return [e_iso], [rev_e_iso]
        e_isos, rev_e_isos = [], []
        for e_H in H.hyper_edges:
            if e_H.id in rev_e_iso:
                continue
            e = hyperedges[i]
            if e_H.label != e.label:
                continue
            compatible = True
            for v in e.get_vertices():
                if v_iso[v] not in e_H.vertices:
                    compatible = False
            if compatible:
                n_e_iso = copy.deepcopy(e_iso)
                n_rev_e_iso = copy.deepcopy(rev_e_iso)
                n_e_iso[e.id] = e_H.id
                n_rev_e_iso[e_H.id] = e.id
                n_e_isos, n_rev_e_isos = self.find_isomorphisms_e(i+1, hyperedges, n_e_iso, n_rev_e_iso, v_iso, rev_v_iso, H)
                for j in range(len(n_e_isos)):
                    e_isos.append(n_e_isos[j])
                    rev_e_isos.append(n_rev_e_isos[j])
        return e_isos, rev_e_isos

    def find_isomorphisms_v(self, i, vertices, v_iso, rev_v_iso, H):
        if i == len(vertices):
            return [v_iso], [rev_v_iso]
        n_v_isos, n_rev_v_isos = [], []
        for v_H in H.vertex_label:
            if v_H in rev_v_iso:
                continue
            v = vertices[i]
            if self.vertex_label[v] == H.vertex_label[v_H]:
                new_v_iso = copy.deepcopy(v_iso)
                new_rev_v_iso = copy.deepcopy(rev_v_iso)
                new_v_iso[v] = v_H
                new_rev_v_iso[v_H] = v
                n_v_isos_temp, n_rev_v_isos_temp = self.find_isomorphisms_v(i + 1, vertices, new_v_iso, new_rev_v_iso, H)
                for j in range(len(n_v_isos_temp)):
                    n_v_isos.append(n_v_isos_temp[j])
                    n_rev_v_isos.append(n_rev_v_isos_temp[j])
        return n_v_isos, n_rev_v_isos

    def add_vertex(self, v_id, v_label):
        if v_id in self.vertex_label and self.vertex_label[v_id] == v_label:
            return
        if v_id in self.vertex_label:
            raise ValueError("Vertex ID exists")
        self.vertex_label[v_id] = v_label
        if v_label in self.vertex_label_bucket:
            self.vertex_label_bucket[v_label] = self.vertex_label_bucket[v_label] + 1
        else:
            self.vertex_label_bucket[v_label] = 1
        self.vertex_count += 1

    def add_hyper_edge(self, vertices, edge_label, weight):
        self.hyper_edges.append(HyperEdge(self.edge_count, vertices, edge_label, weight))
        self.weight += weight
        if edge_label in self.edge_label_bucket:
            self.edge_label_bucket[edge_label] += 1
        else:
            self.edge_label_bucket[edge_label] = 1
        self.edge_count += 1

    def extend_hyper_edge(self, edge_id, vertex_id):
        self.hyper_edges[edge_id].vertices.append(vertex_id)

    def get_new_vertex_id(self):
        return self.vertex_count

    def get_vertex_label(self, id):
        return self.vertex_label[id]

    def get_hyperedge(self, id):
        return self.hyper_edges[id]

    def get_weight(self):
        return self.weight

    def get_hyperedges(self):
        return self.hyper_edges

    def isomorphic(self, h):
        if len(self.hyper_edges) != len(h.hyper_edges):
            return False
        for v in self.vertex_label_bucket:
            if v not in h.vertex_label_bucket:
                return False
            if self.vertex_label_bucket[v] != h.vertex_label_bucket[v]:
                return False
        for edge_label in self.edge_label_bucket:
            if edge_label not in h.edge_label_bucket:
                return False
            if self.edge_label_bucket[edge_label] != h.edge_label_bucket[edge_label]:
                return False
        return self.check_isomorphism(dict(), dict(), h)

    def check_isomorphism(self, isomorphism, rev_isomorphism, h):
        if len(isomorphism) == len(self.vertex_label):#!
            rev_edge_iso = dict()
            for hyper_edge in self.hyper_edges:
                iso_vertices = [isomorphism[v] for v in hyper_edge.vertices]
                new_hyper_edge = HyperEdge(0, iso_vertices, hyper_edge.label)
                find = False
                for iso_hyper_edge in h.hyper_edges:
                    if iso_hyper_edge.id in rev_edge_iso:
                        continue
                    if new_hyper_edge.equals(iso_hyper_edge):
                        rev_edge_iso[iso_hyper_edge.id] = hyper_edge.id
                        find = True
                        break
                if not find:
                    return False
            return True

        for u in self.vertex_label:
            if u not in isomorphism:
                for v in h.vertex_label:
                    if v not in rev_isomorphism and self.vertex_label[u] == h.vertex_label[v]:
                        new_isomorphism = copy.deepcopy(isomorphism)
                        new_rev_isomorphism = copy.deepcopy(rev_isomorphism)
                        new_isomorphism[u] = v
                        new_rev_isomorphism[v] = u
                        if self.check_isomorphism(new_isomorphism, new_rev_isomorphism, h):
                            return True
        return False

    def __str__(self):
        result = ""
        for v in self.vertex_label:
            result += "v " + str(v) + " " + str(self.vertex_label[v]) + "\n"
        for e in self.hyper_edges:
            result += "e " + str(e.vertices) + " " + str(e.label) + " " + str(e.weight) + "\n"
        return result


class DFStuple:
    def __init__(self, vertices, edge_label):
        self.vertices = vertices
        self.edge_label = edge_label

    def add_vertex(self, id, label):
        for v in self.vertices:
            if v >= id:
                return False
        self.vertices[id] = label


class DFScode:

    def __init__(self):
        self.tuples = []
        self.extentions = []
        self.last_edge = -1
        self.last_vertex = -1

    def extend(self, ext):
        self.extentions.append(ext)
        if ext[0] == "a":
            tuple = DFStuple({ext[1]: ext[2]}, ext[3])
            for v in tuple.vertices:
                if v > self.last_vertex:
                    self.last_vertex = v
            self.last_edge += 1
            self.tuples.append(tuple)
        else:
            vertex_id, vertex_label = ext[1], ext[2]
            if vertex_id is None:
                self.last_vertex += 1
                self.tuples[self.last_edge].add_vertex(self.last_vertex, vertex_label)
            else:
                return self.tuples[self.last_edge].add_vertex(vertex_id, vertex_label)

    def get_hypergraph(self):
        h = HyperGraph()
        for t in self.tuples:
            vertices = []
            for v in t.vertices:
                vertices.append(v)
                h.add_vertex(v, t.vertices[v])
            h.add_hyper_edge(vertices, t.edge_label, 0)
        return h

    def get_extensions(self):
        return self.extentions

    def get_last_edge(self):
        return self.last_edge

    def __str__(self):
        result = ""
        for t in self.tuples:
            result += str(t.vertices) + " " + str(t.edge_label) + "\n"
        return result


def load_hypergraphs_from_file(file):
    hypergraphs = []
    hypergraph = None
    for line in file:
        line = line.replace("\n", "")
        if line[0] == 't':
            if hypergraph is not None:
                hypergraphs.append(hypergraph)
            hypergraph = HyperGraph()
        if line[0] == 'v':
            _, v_id, v_label = line.split(" ")
            v_id = int(v_id)
            v_label = int(v_label)
            hypergraph.add_vertex(v_id, v_label)
        if line[0] == 'e':
            temp = line.split(" ")
            vertices = []
            for v in temp[1:len(temp)-2]:
                vertices.append(int(v))
            edge_label = int(temp[len(temp)-2])
            weight = float(temp[len(temp)-1])
            hypergraph.add_hyper_edge(vertices, edge_label, weight)
    if hypergraph is not None:
        hypergraphs.append(hypergraph)
    return hypergraphs


def min_ext(ext1, ext2):
    ext_type1, time1, v_label1, e_label1 = ext1
    ext_type2, time2, v_label2, e_label2 = ext2
    if ext_type1 == "a" and ext_type2 == "e":
        return ext2
    if ext_type1 == "e" and ext_type2 == "a":
        return ext1
    if ext_type1 == "a" and ext_type2 == "a":
        if time1 < time2:
            return ext1
        elif time1 > time2:
            return ext2
        elif e_label1 < e_label2:
            return ext1
        elif e_label1 > e_label2:
            return ext2
        elif e_label1 == e_label2:
            if v_label1 < v_label2:
                return ext1
            else:
                return ext2
    if ext_type1 == "e" and ext_type2 == "e":
        if time1 is None and time2 is None:
            if v_label1 < v_label2:
                return ext1
            else:
                return ext2
        elif time1 is None:
            return ext2
        elif time2 is None:
            return ext1
        elif time1 < time2:
            return ext1
        elif time1 > time2:
            return ext2
        elif time1 == time2:
            if v_label1 < v_label2:
                return ext1
            else:
                return ext2


def get_minimum_ext(extensions):
    result = None
    for ext in extensions:
        if result is None:
            result = ext
        else:
            result = min_ext(result, ext)
    return result


def is_minimum(code):
    h = code.get_hypergraph()
    c = DFScode()
    for i in range(len(code.get_extensions())):
        ext = get_minimum_ext(get_extensions(c, [h]))
        if ext != code.get_extensions()[i]:
            return False
        c.extend(ext)
    return True


def get_extensions(code, hypergraphs):
    extensions = dict()
    progress = 0
    for h in hypergraphs:
        progress += 1
        #print(progress,len(hypergraphs))
        temp_extensions = dict()
        if len(code.tuples) == 0:
            for he in h.hyper_edges:
                for v in he.vertices:
                    if ("a", 0, h.vertex_label[v], he.label) in temp_extensions:
                        w = temp_extensions[("a", 0, h.vertex_label[v], he.label)]
                        temp_extensions[("a", 0, h.vertex_label[v], he.label)] = max(w, he.weight)
                    else:
                        temp_extensions[("a", 0, h.vertex_label[v], he.label)] = he.weight
        else:
            code_graph = code.get_hypergraph()
            v_isos, rev_v_isos, e_isos, rev_e_isos = code_graph.find_isomorphisms(h)
            for i in range(len(v_isos)):
                v_iso, rev_v_iso, e_iso, rev_e_iso = v_isos[i], rev_v_isos[i], e_isos[i], rev_e_isos[i]
                weight = 0
                for e in rev_e_iso:
                    he = h.get_hyperedge(e)
                    weight += he.weight
                he = h.get_hyperedge(e_iso[code.get_last_edge()])
                for v in he.get_vertices():
                    if v not in rev_v_iso:
                        if ("e", None, h.vertex_label[v], None) in temp_extensions:
                            w = temp_extensions[("e", None, h.vertex_label[v], None)]
                            temp_extensions[("e", None, h.vertex_label[v], None)] = max(w, weight)
                        else:
                            temp_extensions[("e", None, h.vertex_label[v], None)] = weight
                    elif rev_v_iso[v] not in code_graph.get_hyperedge(code.get_last_edge()).get_vertices():
                        if ("e", rev_v_iso[v], h.vertex_label[v], None) in temp_extensions:
                            w = temp_extensions[("e", rev_v_iso[v], h.vertex_label[v], None)]
                            temp_extensions[("e", rev_v_iso[v], h.vertex_label[v], None)] = max(w, weight)
                        else:
                            temp_extensions[("e", rev_v_iso[v], h.vertex_label[v], None)] = weight
                for he in h.get_hyperedges():
                    if he.get_id() not in rev_e_iso:
                        for v in he.get_vertices():
                            if v in rev_v_iso:
                                if ("a", rev_v_iso[v], h.vertex_label[v], he.label) in temp_extensions:
                                    w = temp_extensions[("a", rev_v_iso[v], h.vertex_label[v], he.label)]
                                    temp_extensions[("a", rev_v_iso[v], h.vertex_label[v], he.label)] = max(w, weight)
                                else:
                                    temp_extensions[("a", rev_v_iso[v], h.vertex_label[v], he.label)] = weight
        for key in temp_extensions:
            w = temp_extensions[key]
            if key in extensions:
                weight = extensions[key]
                extensions[key] = w * weight
            else:
                extensions[key] = w
    return extensions


def extend_code(code, extensions):
    new_codes = []
    for ext in extensions:
        new_code = copy.deepcopy(code)
        if new_code.extend(ext) is None:
            new_codes.append((new_code,extensions[ext]))
    return new_codes


domains = ["Machine_learning", "Computer_network", "Computer_security", "Data_mining", "Distributed_computing", "Bioinformatics"]
domain = domains[5]
hypergraphs = load_hypergraphs_from_file(open("uncertain_hypergraph_data/" + domain + ".txt").readlines())
out = open("patterns/uncertain/" + domain + ".txt", "w")


found, err, candidates, t = 0, 0, 0, time.time()


def search(code, hypergraphs, min_sup):
    global found, err, t, candidates
    print(time.time()-t)
    print("extending", code.get_extensions())
    extensions = get_extensions(code, hypergraphs)
    print("extended", code.get_extensions())
    codes = extend_code(code, extensions)
    for c in codes:
        candidates += 1
        weight = c[1]
        c = c[0]
        if is_minimum(c) is False:  # Remove this check for getting naive algorithm
            continue
        if weight >= min_sup:
            out.write("t # " + str(found) + "\n")
            out.write(str(c.get_hypergraph())+"\n")
            found += 1
            print("found", c.get_extensions(), weight)
            try:
                search(c, hypergraphs, min_sup)
            except MemoryError:
                err += 1
                pass


delta = 0.0009
min_sup = len(hypergraphs) * delta


t = time.time()
search(DFScode(), hypergraphs, min_sup)
print(delta)
print(err, found, candidates)
print("Runtime:", time.time()-t)
