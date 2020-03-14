#
# Source: https://gist.github.com/joninvski/701720
#


# Step 1: For each node prepare the destination and predecessor
def initialize(graph, source):
    # Stands for destination
    d = {}
    # Stands for predecessor
    p = {}
    for node in graph:
        # We start admiting that the rest of nodes are very very far
        d[node] = float('Inf')
        # For the source we know how to reach
        p[node] = None
    d[source] = 0
    return d, p


def relax(node, neighbour, graph, d, p):
    # If the distance between the node and the neighbour is lower than the one I have now
    if d[neighbour] > d[node] + graph[node][neighbour]:
        # Record this lower distance
        d[neighbour] = d[node] + graph[node][neighbour]
        p[neighbour] = node


def bellman_ford(graph, source):
    d, p = initialize(graph, source)

    # Run this until is converges
    for i in range(len(graph) - 1):
        for u in graph:
            # For each neighbour of u
            for v in graph[u]:
                # Lets relax it
                relax(u, v, graph, d, p)

    # Step 3: check for negative-weight cycles
    for u in graph:
        for v in graph[u]:
            assert d[v] <= d[u] + graph[u][v]

    return d, p
