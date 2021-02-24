# -*- coding: utf-8 -*-


import heapq, itertools


__author__ = 'David Schaller'


class _TreeNode:
    
    def __init__(self, priority, value):
        
        self.priority = priority
        self.value = value
        self.children = []
    
    
    def __eq__(self, o):
        
        return self.priority == o.priority
    
    
    def __le__(self, o):
        
        return self.priority <= o.priority
    
    
    def __lt__(self, o):
        
        return self.priority < o.priority


def karmarkar_karp(iterable):
    
    heap = []
    
    for item in iterable:
        if isinstance(item, int):
            heapq.heappush(heap, (_TreeNode(-item, item)))
        else:
            heapq.heappush(heap, (_TreeNode(-len(item), item)))
    
    while len(heap) > 1:
        x = heapq.heappop(heap)
        y = heapq.heappop(heap)
        x.children.append(y)
        x.priority -= y.priority
        heapq.heappush(heap, x)
    
    root = heap[0]
    
    # construct the actual partition via 2-coloring of the tree
    part = [[],[]]
    stack = [(root, 0)]
    
    while stack:
        node, color = stack.pop()
        part[color].append(node.value)
        
        for child in node.children:
            stack.append( (child, (color+1)%2) )
    
    return part


def balanced_coarse_graining(partition):
    
    return [[item for item in itertools.chain.from_iterable(set_i)]
            for set_i in karmarkar_karp(partition)]


if __name__ == '__main__':
    
    numbers = [8, 7, 6, 5, 4]
    print(karmarkar_karp(numbers))
    
    part = [[3,2,1], [4,5,6,7], [8,9], [10,11,12,13,14], [15]]
    print(balanced_coarse_graining(part))