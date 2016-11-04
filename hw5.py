import custom_utils as utils

def greedy_clique_partition(adj_list):
    i = 0
    Q = set()
    clusters = []
    V = set()
    
    for node in adj_list:
        V.add(node)
    
    while len(V-Q) != 0:
        i += 1
        C = set()
        V_prime = V - Q
        while len(V_prime) > 0:
            m = get_max_degree(V_prime, adj_list)
            C.add(m)
            N = set()
            for neighbor in adj_list[m]:
                N.add(neighbor)
            V_prime = V_prime & N
        Q = Q | C
        clusters.append(C)
        clusters.sort(key=lambda x: len(x), reverse=True)
    return clusters

def get_clique_dict(adj_ls):
    #returns a dictionary with pairs of form (node: clique number)
    d = {}
    c = greedy_clique_partition(adj_ls)
    for i in range(len(c)):
        for node in c[i]:
            d[node] = i
    return d

    
def get_max_degree(S, adj_list):
    #helper function for greedy_clique_partition()
    #finds the node with max degree of the set S
    return max(S, key = lambda n: len(adj_list[n]))
    
def main():
    example_graph,nodes,edges = utils.parse_input('example.txt','\t')
    c_d = get_clique_dict(example_graph.get_adj_ls())
    nodeAttrs = example_graph.GSnodeAttrs
    bg = 'background_color'
    
    for n in nodeAttrs:
        if c_d[n] == 0:
            nodeAttrs[n][bg] = 'red'
        elif c_d[n] == 1:
            nodeAttrs[n][bg] = 'blue'
        elif c_d[n] == 2:
            nodeAttrs[n][bg] = 'green'
        elif c_d[n] == 3:
            nodeAttrs[n][bg] = 'cyan'
        elif c_d[n] == 4:
            nodeAttrs[n][bg] = 'magenta'

    example_graph.GSnodeAttrs = nodeAttrs

    example_graph.uploadGraph()
    
    
if __name__ == '__main__':
    main()
