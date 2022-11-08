import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

import Base: coalesce
using DataFrames
using GLMakie, GraphMakie, NetworkLayout
using Graphs
import StatsBase: countmap, sample

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
  while i!=j
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
function A(P::SimpleDiGraph, a=fill(0, nv(P)), maxiter=100000)
  leaves = filter(n -> isleaf(P, n), vertices(P))
  ready = falses(size(a))
  a[leaves] .= 0
  ready[leaves] .= true

  _bothready(v) = bothready(P, ready, v)
  i = findfirst(_bothready, eachindex(ready))
  iter = 0
  while !isnothing(i) && iter<maxiter
    c = children(P, i)
    a[i] += a[c[1]] + a[c[2]] + 2
    ready[i] = true
    ready[c] .= false
    i = findfirst(_bothready, eachindex(ready))
    iter += 1
  end
  iter>=maxiter && @warn("maxiter reached.")

  return a.+1
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

function moran(n, T, mu=0.0, earlystop=false)
  t = 1
  moves = Vector{Tuple{Int,Int}}(undef, T)
  v = collect(1:n)
  g = n
  while t<=T 
    i,j = sample(eachindex(v), 2, replace=false)
    if rand() < mu
      g += 1
      v[j] = g
    else
      v[j] = v[i]
    end
    moves[t] = (i,j)
    t += 1
  end

  return v, moves
end

function coalesce(P::SimpleDiGraph, i, j)
  add_vertex!(P)
  v = nv(P)
  i != -1 && add_edge!(P, i, v)
  add_edge!(P, j, v)
  return v
end

function find_ancestry(n, moves)
  P = SimpleDiGraph(n)
  current_gen = collect(1:n)
  active = n
  for (i,j) in Iterators.reverse(moves)
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

  return P, current_gen
end

function yule(n)
  P = SimpleDiGraph(1)
  leaves = [1]
  k = 1
  while k<n
    i = rand(eachindex(leaves))
    v = leaves[i]
    add_vertex!(P)
    add_vertex!(P)
    v1, v2 = last(vertices(P), 2)
    add_edge!(P, v1, v)
    add_edge!(P, v2, v)
    deleteat!(leaves, i)
    append!(leaves, [v1,v2])
    k += 1
  end
  P
end


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

function birthdeath(n, T, d, b=1.0; N=0)
  P = SimpleDiGraph(n)
  t = 1
  # moves = Vector{Tuple{Int,Int}}(undef, T)
  leaves = collect(vertices(P))
  N>0 && sizehint!(leaves, N)
  while t<=T && (N==0 || length(leaves)<N)
    isempty(leaves) && break 
    i = rand(eachindex(leaves))
    if rand()<b
      # insert!(pop, i+1, pop[i])
      add_vertex!(P)
      add_vertex!(P)
      v1,v2 = last(vertices(P), 2)
      add_edge!(P, v1, leaves[i])
      add_edge!(P, v2, leaves[i])
      insert!(leaves, i+1, v1)
      insert!(leaves, i+2, v2)
      deleteat!(leaves, i)
      # moves[t] = (i,i+1)
    end
    j = rand(eachindex(leaves))
    if rand()<d # die
      deletevertex!(P, leaves, j)
    end
    t += 1
  end

  root = findfirst(v->isnothing(parent(P,v)), vertices(P))

  M = adjacency_matrix(P)
  tmp = M[1, :]
  M[1, :] .= M[root, :]
  M[root, :] .= tmp
  tmp = M[:, 1]
  M[:, 1] .= M[:, root]
  M[:, root] .= tmp

  return SimpleDiGraph(M)
end

function moran2(n, T)
  P = SimpleDiGraph(n)
  t = 1
  # moves = Vector{Tuple{Int,Int}}(undef, T)
  pop = collect(1:n)
  vert = collect(vertices(P))
  while t<=T
    isempty(vert) && break 
    i = rand(eachindex(vert))
    # insert!(pop, i+1, pop[i])
    add_vertex!(P)
    add_vertex!(P)
    v1,v2 = last(vertices(P), 2)
    add_edge!(P, v1, vert[i])
    add_edge!(P, v2, vert[i])
    insert!(vert, i+1, v1)
    insert!(vert, i+2, v2)
    deleteat!(vert, i)
    # moves[t] = (i,i+1)

    deletevertex!(P, vert, rand(eachindex(vert)))
    t += 1
  end

  root = findfirst(v->isnothing(parent(P,v)), vertices(P))

  M = adjacency_matrix(P)
  tmp = M[1, :]
  M[1, :] .= M[root, :]
  M[root, :] .= tmp
  tmp = M[:, 1]
  M[:, 1] .= M[:, root]
  M[:, root] .= tmp

  return SimpleDiGraph(M)
end

function maximally_unbalanced(height)
  P = SimpleDiGraph(1)
  h = 0
  attach_to = 1
  while h<height
    add_vertex!(P)
    add_vertex!(P)
    v1,v2 = last(vertices(P), 2)
    add_edge!(P, v1, attach_to)
    add_edge!(P, v2, attach_to)
    attach_to = v1
    h += 1
  end
  P
end    

function maximally_balanced(height)
  P = SimpleDiGraph(1)
  h = 1
  attach_to = [1]
  while h<height
    v = popfirst!(attach_to)
    add_vertex!(P)
    add_vertex!(P)
    v1,v2 = last(vertices(P), 2)
    add_edge!(P, v1, v)
    add_edge!(P, v2, v)
    append!(attach_to, [v1, v2])
    h = log2(nv(P)+1)
  end
  P
end
      
