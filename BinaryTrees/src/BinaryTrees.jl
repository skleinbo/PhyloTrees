module BinaryTrees

export adjacency_matrix, BinaryTree, child!, left!, right!, sibling, isleftchild, isrightchild

import AbstractTrees
import AbstractTrees: children, nodevalue, parent, ParentLinks, PreOrderDFS, StoredParents
import Base: show
import SparseArrays: spzeros

mutable struct BinaryTree{T}
    val::T
    parent::Union{Nothing, BinaryTree{T}}
    left::Union{Nothing, BinaryTree{T}}
    right::Union{Nothing, BinaryTree{T}}
    function BinaryTree{T}(v=zero(T), parent=nothing, left=nothing, right=nothing) where T
        new{T}(v, parent, left, right)
    end
end
BinaryTree(v::T) where T = BinaryTree{T}(v, nothing, nothing, nothing)

function AbstractTrees.children(t::BinaryTree) 
    if isnothing(t.left) && isnothing(t.right)
        return ()
    elseif isnothing(t.right)
        return (t.left,)
    elseif isnothing(t.left)
        return (t.right,)
    else
        return (t.left, t.right)
    end
end
AbstractTrees.nodevalue(t::BinaryTree) = t.val
AbstractTrees.parent(t::BinaryTree) = t.parent
AbstractTrees.ParentLinks(::Type{<:BinaryTree}) = AbstractTrees.StoredParents()

function child!(t::BinaryTree{T}, v::T, location::Symbol) where T
    newnode = BinaryTree(v)
    newnode.parent = t
    if location===:left
        t.left = newnode
    elseif location===:right
        t.right = newnode
    else
        throw(ArgumentError("location must be one of :left, :right"))
    end
    return newnode
end

left!(t::BinaryTree{T}, v::T) where T = child!(t, v, :left)
right!(t::BinaryTree{T}, v::T) where T = child!(t, v, :right)

function sibling(t::BinaryTree)
    p = parent(t)
    isnothing(p) && return nothing
    c = children(p)
    length(c) == 2 || return nothing
    return c[1]===t ? c[2] : c[1]
end

isleftchild(t::BinaryTree)  = !isnothing(parent(t)) && t===parent(t).left
isrightchild(t::BinaryTree) = !isnothing(parent(t)) && t===parent(t).right

function show(io::IO, ::MIME"text/plain", t::BinaryTree{T}) where T
    print(io, "BinaryTree{$T} $(nodevalue(t)) $(objectid(t)) with $(length(children(t))) children")
    if isnothing(parent(t))
        print(io, " and no parent.")
    else
        print(io, ".")
    end
    print(io, "\n")
    nothing
end

function adjacency_matrix(t::BinaryTree)
    n = 0
    for _ in PreOrderDFS(t)
        n += 1
    end

    A = spzeros(Int, n,n)
    id = IdDict(t=>1)
    i = 2
    for v in PreOrderDFS(t)
        p = parent(v)
        isnothing(p) && continue
        if !(v in keys(id))
            id[v] = i
            A[id[p], i] = 1
            i += 1
        end
    end

    return A
end

end # MODULE