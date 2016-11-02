
def greedy_clique_partition(adj_list):
    i = 0
    Q = set()
    clusters = []
    V = set()
    
    for node in adj_list:
        V.add(node)
        
    while len(V-Q) == 0:
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
    return clusters


    
def get_max_degree(S, adj_list):
    return max(S, key = lambda n: len(adj_list[n]))
    
    
    

        