# calculate sub-tree height...

# define function to test if all array elements equal
allequal(x) = all(y->y==x[1],x)

# should be called treeheight or something really
function mrca_fullstore(ancestry::Array{Int64,2}, leaves::Array{Int64,1})
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
