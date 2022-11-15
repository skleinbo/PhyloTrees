import AbstractTrees: children, Leaves, nodevalue, parent, print_tree, PreOrderDFS, PostOrderDFS
import Base: coalesce
using DataFrames
using StatsBase

include("BinaryTrees.jl")
using .BinaryTrees

"Distance between a node `i` and its ancestor `j` on a tree."
function dist(i, j)
    d = 0
    while i != j
        i = parent(i)
        isnothing(i) && throw(ErrorException("$i is not a descendent of $j."))
        d += 1
    end
    return d
end

## Phylogenetic Observables

function empty!(P::BinaryTree{Vector{T}}) where {T}
    for v::BinaryTree{Vector{T}} in PreOrderDFS(P)
        v.val::Vector{T} .= zero(T)
    end
    nothing
end

function A!(P::BinaryTree{Vector{Int}})
    for v::BinaryTree{Vector{Int}} in PreOrderDFS(P)
        v.val[1] += 1
        p = parent(v)
        while !isnothing(p)
            p.val[1] += 1
            p = parent(p)
        end
    end
    return P
end

function C!(P::BinaryTree{Vector{Int}})
    empty!(P)
    A!(P)
    for v::BinaryTree{Vector{Int}} in PreOrderDFS(P)
        v.val[2] += v.val[1]
        p = parent(v)
        while !isnothing(p)
            p.val[2] += v.val[1]
            p = parent(p)
        end
    end
    return P
end

function treevalues!(P::BinaryTree{Vector{T}}) where {T}
    C!(P)
    V = map(x -> x.val, PreOrderDFS(P))
    return ntuple(i -> getindex.(V, i), length(V[1]))
end

"""
Simulate the coalescent process for n genes.

Return a directed graph with edges pointing towards the root.
"""
function coalescent(n)
    P = [BinaryTree([0, 0]) for _ in 1:n] # start with n vertices; store A,C
    while n > 1
        i, j = sample(1:n, 2, replace=false, ordered=true) #sample two distinct vertices to coalesce
        @inbounds l, r = P[i], P[j]
        v = BinaryTree([0, 0]) # newly added parental node
        v.left = l # connect sampled nodes to parental node
        v.right = r
        l.parent = v
        r.parent = v

        # replace nodes i,j in active nodes list with
        # new node and current last node.
        # shorten list of admissible draws by one.
        P[i] = v
        P[j] = P[n]
        n -= 1
    end
    return P[1]
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
  Construct a Yule tree with `n` tips.

  At each time step a bifurcation event happens at any of the
  current tips with equal probability.
"""
function yule(n)
    P = [BinaryTree([0, 0])]
    k = 1
    while length(P) - k + 1 < n
        i = rand(k:lastindex(P))
        v = P[i]
        l = left!(v, [0, 0])
        r = right!(v, [0, 0])
        append!(P, (l, r))
        P[i] = P[k]

        k += 1
    end
    P[1]
end

"""
  birthdeath(n, T, d, b=1.0; N=0)

  Start with a population of size `n`.
  At each timestep an individual duplicates with 
  probability `b`, and one dies with probability `d`.

  Run for `T` timesteps if `N=0`, or until population has
  reached size `N`; whatever happens first.

  Lineages that die out are pruned automatically.

  Return an ancestral tree, or a vector of trees if the process hasn't coalesced.
"""
function birthdeath(n, T, d, b=1.0; N=0)
    P = [BinaryTree([0, 0]) for _ in 1:n]
    t = 1
    N > 0 && sizehint!(P, N)
    k = 1
    while t <= T && (N == 0 || length(P) - k + 1 < N)
        k > length(P) && break
        if rand() < b
            i = rand(k:lastindex(P))
            v = P[i]
            l = left!(v, [0, 0])
            r = right!(v, [0, 0])
            push!(P, r)
            P[i] = l
        end
        if rand() < d # die
            j = rand(k:lastindex(P))
            w = P[j]
            P[j] = P[k]
            k += 1

            p = parent(w)
            if !isnothing(p)
                pp = parent(p)
                sib = sibling(w)
                if !isnothing(pp)
                    if pp.left === p
                        pp.left = sib
                    else
                        pp.right = sib
                    end
                    sib.parent = pp
                else
                    sib.parent = nothing
                end
            end
        end
        t += 1
    end

    return union(AbstractTrees.getroot.(P[k:end]))
end

"""
  Directly construct an ancestral tree for a Moran process of
    `n` individuals.

  Lineages that die out are pruned automatically.

  See also: @ref(`birth_death`)
"""
moran2(n, T) = birthdeath(n, T, 1.0)

"""
  Return a fully imbalanced binary tree of given height.
"""
function maximally_unbalanced(height)
    P = BinaryTree([0, 0])
    h = 0
    attach_to = P
    while h < height
        right!(attach_to, [0, 0])
        attach_to = left!(attach_to, [0,0])
        h += 1
    end

    return P
end

"""
  Return a fully balanced binary tree of given height.
"""
function maximally_balanced(height)
    P = [BinaryTree([0, 0])]
    root = P[1]
    h = 1
    n = 1
    while h < height
        attach_to = popfirst!(P)
        push!(P, left!(attach_to, [0, 0]))
        push!(P, right!(attach_to, [0, 0]))

        n += 2
        h = log2(n + 1)
    end
    return root
end

"""
  Combine vectors with values for A and C into a dataframe, averaging C/A for every A.
"""
function to_mean_dataframe(A, C)
    df = DataFrame(; A, C)
    combine(groupby(df, :A), :A => mean => :a, [:A, :C] => ((a, c) -> mean(c ./ a)) => :covera)
end

"""
Simulate a "fluctating coalescent" process for n genes.

Return a directed graph with edges pointing towards the root.
"""
function fluctuating_coalescent(n, w=randn(n) .^ 2; default_value=[0, 0])
    P = [BinaryTree(copy(default_value)) for _ in 1:n]
    ij = [0,0]
    while n  > 1
        # i, j = wsample(1:n, @view(w[1:n]), 2, replace=false, ordered=true) #sample two distinct vertices to coalesce
        i,j = StatsBase.sample!(1:n, Weights(view(w, 1:n)), ij)
        l, r = P[i], P[j]
        v = BinaryTree(copy(default_value))
        v.left = l
        v.right = r
        l.parent = v
        r.parent = v

        P[i] = v
        P[j] = P[n]
        w[i] = max(w[i], w[j])
        n -= 1
    end

    return P[1]
end
