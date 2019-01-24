import AbstractTrees: children, printnode, print_tree, Leaves

# define data structure for trees:
struct TreeNode
    generation::Int64
    index::Int64
    children::Array{TreeNode3,1}
end

# make a tree:
node21 = TreeNode(2,1, TreeNode[])
node11 = TreeNode(1,1, [node21])

# method of AbstractTrees::children for type TreeNode:
children(node::TreeNode) = something(node.children, [])
children(node11)

# method of AbstractTrees::printnode for type TreeNode, prints gen.index:
printnode(io::IO, node::TreeNode) = print(io, node.generation, '.', node.index)
print_tree(node11)
