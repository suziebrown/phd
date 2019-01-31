# calculate sub-tree height...

# define function to test if all array elements equal
allequal(x) = all(y->y==x[1],x)

# should be called treeheight or something really
function mrca_fullstore(ancestry::Array{UInt16,2}, leaves::Array{UInt16,1})
    T = size(ancestry)[1]
    mrca = T
    for t in 1:T
        leaves = ancestry[T-t+1, leaves]
        if allequal(leaves)
            mrca = t
            break
        end
    end
    return mrca
end

#= for tree-stored ancestry version, just need to change tht type of the
argument, and use the 'parent' property of the selected nodes in place of
the matrix indexing =#
