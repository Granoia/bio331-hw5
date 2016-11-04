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

    #example_graph.uploadGraph()



    #Cross linking: this technique works by covalently bonding proximal proteins together by use of a chemical agent. Pretty straight forward. If two proteins are interacting they get stuck together.

    CL_graph, clnodes, cledges = utils.parse_input('MI_0030_cross-linking_study.txt','\t',isWeighted=True)
    CL_cliques = greedy_clique_partition(CL_graph.better_adj_ls())
    print(CL_cliques[0])
    #largest clique contains:
    #Transforming growth factor beta receptor type 3
    #    most significant p value GO term (.396): epicardium-derived cardiac fibroblast cell development
    #          i.e. this is non-significantly implicated in maturation of cells that comprise some type of connective tissue in the heart
    #TGF-beta receptor type 2
    #    most significant p value GO term (.396): positive regulation of tolerance induction to self antigen
    #        the description is sparse but I imagine this term is for things that play a part in keeping immune cells from targetting the body
    #Serine/threonine-protein kinase receptor R3
    #    most significant p value GO term (1.00): blood vessel endothelial cell proliferation involved in sprouting angiogenesis
    #        this is the first term on the list but all of them are p = 1.00 so that's a little mysterious. I guess I can talk a little about the function of the term itself. Angiogenesis is creation of blood vessels. Blood vessels have endothelial cells. When you create blood vessels you gotta proliferate some endothelial cells.



    
    #Nuclear Magnetic Resonance (NMR): I would probably need to know more about wave physics to give a good explanation but the gist of it is,
    #you expose your compound to a powerful magnetic field, and certain nuclei (the standard measures tends to be for H-1 or C-13) in the compound will resonate according to the energy of the field.
    #you can detect this resonance to get a signature that can tell you things about the structure of the molecule.
    #fancy versions of NMR such as NOESY (Nuclear Overhauser Effect Spectroscopy) can even tell you when certain nuclei in certain amino acids are within ~5 angstroms of each other, which can give you
    #information about how proteins are interacting.
    NMR_graph, nmrnodes, nmredges = utils.parse_input('MI_0077_nuclear_magnetic_resonance.txt','\t',isWeighted=True)
    NMR_cliques = greedy_clique_partition(NMR_graph.better_adj_ls())
    print(NMR_cliques[0])
    #largest clique contains:
    #Cytoplasmic protein NCK2
    #    most significant p value GO term (1.00): negative regulation of transcription from RNA polymerase II promoter in response to endoplasmic reticulum stress
    #        Again, this is the first one on the list but they all have p = 1.00 so I'm not really sure what to make of that. I guess that's a sign that the protein is poorly characterized which makes sense given that it's called 'Cytoplasmic protein NCK2'. Negative regulation of transcription from RNA polymerase II in response to ER stress would make sense, since ER stress is bound to cause changes to the transcriptome.
    #Polyubiquitin-C
    #    most significant p value GO term (.792): ubiquitin homeostasis
    #        non-significantly involved in maintinence of ubiquitin I guess. Ubiquitin is a protein that tags other proteins for degredation.



    

    #X-ray crystallography: The biochemical gold standard for protein characterization.
    #First you figure out the conditions under which your protein crystallizes (this can take awhile) and then you shoot x-rays at the crystal. This produces a diffraction pattern. Feed the pattern into
    #a computer program and you can end up with a beautifully accurate characterization of your protein's structure! The detail of this technique is fine enough that not only can you see /that/ two proteins
    #are interacting, you can develop a robust picture of /how/. 
    XRC_graph, xrcnodes, xrcedges = utils.parse_input('MI_0114_x-ray_crystallography.txt','\t',isWeighted=True)
    XRC_cliques = greedy_clique_partition(XRC_graph.better_adj_ls())
    print(XRC_cliques[0])
    #largest clique contains:
    #Cellular tumor antigen p53
    #   most significant p value GO term (.729): oligodendrocyte apoptotic process
    #       an oligodendrocyte is a type of cell that forms the myelin sheath for neurons. apoptosis is programmed cell death.
    #Ubiquitin carboxyl-terminal hydrolase 7
    #   most significant p value GO term (1.00): maintinence of DNA methylation
    #       DNA methylation appears to somehow regulation transcription but I don't think there's a consensus on exactly how. I think in a genomics class I've been told that methylation usually causes lower rates of transcription.
    #E3 ubiquitin-protein ligase Mdm2
    #   most significant p value GO term (.396): cellular response to vitamin B1
    #       cells need to respond to chemotaxis in order to effectively regulate their metabolism and other aspects of homeostasis. This is presumably one piece of this process.

    
    
if __name__ == '__main__':
    main()
