# Function for generating list of network edges (to, from) from specified no.
# of random graphs

    # Args
        # x = number of graphs to generate
        # n = number of nodes
        # pm = probability or number of edges
        # type = probability or number (default is "gnm")
        # directed = are graphs directed? (default = FALSE)
        # loops = can there be self loops? (default = FALSE)


gen_rand_edges <- function(x, n, pm, type="gnm", directed=FALSE, loops=FALSE){
    rand_edges <- list()
    for(i in 1:x){
        a <- erdos.renyi.game(n=n, p.or.m = pm, type=type, directed=directed,
                              loops=loops)
        rand_edges[[i]]<- get.edgelist(a)
        rm(a)
    }
    return(rand_edges)
}
