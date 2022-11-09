import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

import Base: coalesce
using DataFrames
using GLMakie, GraphMakie, NetworkLayout
using Graphs
using StatsBase

children(P::SimpleDiGraph, v) = inneighbors(P, v)
nchildren(P::SimpleDiGraph, v) = length(children(P, v))
has_children(P, v) = nchildren(P, v) > 0
isleaf(P, v) = !has_children(P, v)

"Return the parent vertex of `v` in a tree `P`, or `nothing` if `v` is the root"
function parent(P::SimpleDiGraph, v)
    on = outneighbors(P, v)
    if isempty(on)
        return nothing
    end
    return on[1]
end

"Distance between a node `i` and its parent `j` on tree `P`."
function dist(P, i, j)
    d = 0
    while i != j
        i = parent(P, i)
        isnothing(i) && throw(ErrorException("$i is not a descendent of $j."))
        d += 1
    end
    return d
end

"Helper function for `A`"
function bothready(P, ready, v)
    c = children(P, v)
    return !isempty(c) && ready[c[1]] && ready[c[2]]
end

## Phylogenetic Observables
## TODO: Make this simpler
function A(P::SimpleDiGraph, a=fill(0, nv(P)), maxiter=100000)
    leaves = filter(n -> isleaf(P, n), vertices(P))
    ready = falses(size(a))
    a[leaves] .= 0
    ready[leaves] .= true

    _bothready(v) = bothready(P, ready, v)
    i = findfirst(_bothready, eachindex(ready))
    iter = 0
    while !isnothing(i) && iter < maxiter
        c = children(P, i)
        a[i] += a[c[1]] + a[c[2]] + 2
        ready[i] = true
        ready[c] .= false
        i = findfirst(_bothready, eachindex(ready))
        iter += 1
    end
    iter >= maxiter && @warn("maxiter reached.")

    return a .+ 1
end

function C(P::SimpleDiGraph, a=A(P))
    c = copy(a)
    for v in vertices(P)
        children = neighborhood(P, v, nv(P), dir=:in) # find all children
        popfirst!(children) # remove v itself
        isempty(children) && continue
        # dists = dist.(Ref(P), children, v)
        @views c[v] = sum(a[children]) + a[v]
    end
    return c
end

"""
Simulate the coalescent process for n genes.

Return a directed graph with edges pointing towards the root.
"""
function coalescent(n)
    P = SimpleDiGraph(n)
    V = collect(1:n) # start with n vertices
    while length(V) > 1
        i, j = sample(eachindex(V), 2, replace=false) #sample two distinct vertices to coalesce
        add_vertex!(P) # introduce ancestral node
        v = nv(P) # index of newly added node
        add_edge!(P, V[i], v) # connect samples nodes to parental node
        add_edge!(P, V[j], v)
        deleteat!(V, sort([i, j])) # remove coalesced nodes from list of "active" node

        push!(V, v) # mark new node as active
    end
    P
end

## Forwards Moran process ##

"""
  Simulate a forward Moran process.

  Start with a population 1...n.

  Return a population (i_1, i_2,...,i_n), and
  a list of moves (i,j) indicating individual at i 
  replaced the one at j.
"""
function moran(n, T)
    t = 1
    moves = Vector{Tuple{Int,Int}}(undef, T)
    v = collect(1:n)
    while t <= T
        i, j = sample(eachindex(v), 2, replace=false)
        v[j] = v[i]
        moves[t] = (i, j)
        t += 1
    end

    return v, moves
end

"""
  Helper function to `find_ancestry`.
  
  Introduce a new vertex, and connect i and j to it.
  Unless i=-1 in which case only j is connected.
"""
function coalesce(P::SimpleDiGraph, i, j)
    add_vertex!(P)
    v = nv(P)
    i != -1 && add_edge!(P, i, v)
    add_edge!(P, j, v)
    return v
end

"""
  Construct an ancestral tree from a list of moves.

  Lineages that die out are automatically pruned.

  Return a tree with the MRCAs as roots.
"""
function find_ancestry(n, moves)
    P = SimpleDiGraph(n)
    current_gen = collect(1:n)
    active = n
    for (i, j) in Iterators.reverse(moves)
        active == 1 && break
        current_gen[j] == -1 && continue
        if current_gen[i] == -1
            active += 1
            current_gen[i] = current_gen[j]
        else
            newv = coalesce(P, current_gen[i], current_gen[j])
            current_gen[i] = newv
        end
        current_gen[j] = -1
        active -= 1
    end

    return P
end

"""
  Construct a Yule tree with `n` tips.

  At each time step a bifurcation event happens at any of the
  current tips with equal probability.
"""
function yule(n)
    P = SimpleDiGraph(1)
    leaves = [1]
    k = 1
    while k < n
        i = rand(eachindex(leaves))
        v = leaves[i]
        add_vertex!(P)
        add_vertex!(P)
        v1, v2 = last(vertices(P), 2)
        add_edge!(P, v1, v)
        add_edge!(P, v2, v)
        deleteat!(leaves, i)
        append!(leaves, [v1, v2])
        k += 1
    end
    P
end

"""
  deletevertex!(P, vert, i)

  Helper function for `birth_death`.
  
  Deletes vertex `vert[i]` and its parent, and connects the now
  lose subtree back to the tree.
 
  # Explanation:
  Removing a vertex from an ancestral tree renders
  its parental node superflous. We remove it as well, taking
  care to connect the other child node to the grandparent node if
  it exists.

  `vert` is a list of references to some (all) vertices in `P`. When removing a 
  vertex from `P`, vertices get renumbered and the references become 
  invalid. This is corrected by keeping track of the renumbering and updating `vert`
  accordingly.
"""
function deletevertex!(P, vert, i)
    v = vert[i]
    p = parent(P, v)
    pp = nothing
    if !isnothing(p)
        pp = parent(P, p)
        if !isnothing(pp)
            c = children(P, p)
            c = c[1] == v ? c[2] : c[1]
            add_edge!(P, c, pp)
        end
        lastv = nv(P)
        ilastv = findfirst(==(lastv), vert)
        if !isnothing(ilastv)
            vert[ilastv] = p
        end
        if v == lastv
            v = p
            i = ilastv
        end
        rem_vertex!(P, p)
    end
    lastv = nv(P)
    ilastv = findfirst(==(lastv), vert)
    if !isnothing(ilastv)
        vert[ilastv] = v
    end
    rem_vertex!(P, v)
    deleteat!(vert, i)
end


"""
  birthdeath(n, T, d, b=1.0; N=0)

  Start with a population of size `n`.
  At each timestep an individual duplicates with 
  probability `b`, and one dies with probability `d`.

  Run for `T` timesteps if `N=0`, or until population has
  reached size `N`; whatever happens first.

  Lineages that die out are pruned automatically.

  Return an ancestral tree.
"""
function birthdeath(n, T, d, b=1.0; N=0)
    P = SimpleDiGraph(n)
    t = 1
    leaves = collect(vertices(P))
    N > 0 && sizehint!(leaves, N)
    while t <= T && (N == 0 || length(leaves) < N)
        isempty(leaves) && break
        i = rand(eachindex(leaves))
        if rand() < b
            add_vertex!(P)
            add_vertex!(P)
            v1, v2 = last(vertices(P), 2)
            add_edge!(P, v1, leaves[i])
            add_edge!(P, v2, leaves[i])
            insert!(leaves, i + 1, v1)
            insert!(leaves, i + 2, v2)
            deleteat!(leaves, i)
        end
        j = rand(eachindex(leaves))
        if rand() < d # die
            deletevertex!(P, leaves, j)
        end
        t += 1
    end

    ## Because the root is generally not vertex no. 1, shuffle it there.
    ## Needed in conjunction with GraphMakie and Buchheim tree layout.
    root = findfirst(v -> isnothing(parent(P, v)), vertices(P))

    M = adjacency_matrix(P)
    tmp = M[1, :]
    M[1, :] .= M[root, :]
    M[root, :] .= tmp
    tmp = M[:, 1]
    M[:, 1] .= M[:, root]
    M[:, root] .= tmp

    return SimpleDiGraph(M)
end

"""
  Directly construct an ancestral tree for a Moran process of
    `n` individuals.

  Lineages that die out are pruned automatically.

  See also: @ref(`birth_death`)
"""
function moran2(n, T)
    P = SimpleDiGraph(n)
    t = 1
    # moves = Vector{Tuple{Int,Int}}(undef, T)
    pop = collect(1:n)
    vert = collect(vertices(P))
    while t <= T
        isempty(vert) && break
        i = rand(eachindex(vert))
        # insert!(pop, i+1, pop[i])
        add_vertex!(P)
        add_vertex!(P)
        v1, v2 = last(vertices(P), 2)
        add_edge!(P, v1, vert[i])
        add_edge!(P, v2, vert[i])
        insert!(vert, i + 1, v1)
        insert!(vert, i + 2, v2)
        deleteat!(vert, i)
        # moves[t] = (i,i+1)

        deletevertex!(P, vert, rand(eachindex(vert)))
        t += 1
    end

    root = findfirst(v -> isnothing(parent(P, v)), vertices(P))

    M = adjacency_matrix(P)
    tmp = M[1, :]
    M[1, :] .= M[root, :]
    M[root, :] .= tmp
    tmp = M[:, 1]
    M[:, 1] .= M[:, root]
    M[:, root] .= tmp

    return SimpleDiGraph(M)
end

"""
  Return a fully imbalanced binary tree of given height.
"""
function maximally_unbalanced(height)
    P = SimpleDiGraph(1)
    h = 0
    attach_to = 1
    while h < height
        add_vertex!(P)
        add_vertex!(P)
        v1, v2 = last(vertices(P), 2)
        add_edge!(P, v1, attach_to)
        add_edge!(P, v2, attach_to)
        attach_to = v1
        h += 1
    end
    P
end

"""
  Return a fully balanced binary tree of given height.
"""
function maximally_balanced(height)
    P = SimpleDiGraph(1)
    h = 1
    attach_to = [1]
    while h < height
        v = popfirst!(attach_to)
        add_vertex!(P)
        add_vertex!(P)
        v1, v2 = last(vertices(P), 2)
        add_edge!(P, v1, v)
        add_edge!(P, v2, v)
        append!(attach_to, [v1, v2])
        h = log2(nv(P) + 1)
    end
    P
end

"""
  Combine vectors with values for A and C into a dataframe, averaging C/A for every A.
"""
function to_mean_dataframe(A, C)
    df = DataFrame(; A, C)
    combine(groupby(df, :A), :A => mean => :a, [:A, :C] => ((a, c) -> mean(c ./ a)) => :covera)
end
