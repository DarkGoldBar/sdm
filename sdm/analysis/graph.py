# -*- coding: utf-8 -*-
def BFSTraversal(adj, start: int, data: list = None) -> list:
    """以图的一个点为起点，对这个点的每个分支做并行广度优先搜索遍历，用于判断基团是否拓扑相同
    集合overlap 用于处理遇到环的情况，同层遇到同一节点则同时抛弃这个节点

    Args:
        adj: (dict-like) Graph.adj 临接字典
        start: 起点的key
        data: (list | None) 每个节点的数据，如果无数据则返回节点id

    Return:
        [[(0,), (1, 6), (2, 17, 18, 19), (7, 9), (8, 20), (14, 15, 16)],
         [(4,), (3, 12), (2, 11)],
         [(13,)]]
    """
    data = list(adj) if data is None else data
    history = [[] for x in adj[start]]
    visited = {start}
    current = [{x} for x in adj[start]]
    while any(current):
        for h, c in zip(history, current):
            if c:
                h.append(tuple(data[i] for i in sorted(c)))
        new_cur_all = [i for c in current for i in c]
        visited |= set(new_cur_all)
        overlap = set(i for i in new_cur_all if new_cur_all.count(i) > 1)
        current = [c - overlap for c in current]
        current = [set(ii for i in c for ii in adj[i]) - visited for c in current]
    return history


def isRotatableBond(adj, i: int, j: int, at: list) -> bool:
    adj_i = adj[i]
    adj_j = adj[j]  # noqa
    # 末端原子、全氢、全氟
    if (len(adj_i) == 1) or \
        all(at[x] == 'H' for x in adj_i if x != j) or \
        all(at[x] == 'F' for x in adj_i if x != j):
        return False
    # 羟基内部的要算
    # O_j---H
    # j_is_oh = False
    # if at[j] == 'O' and len(adj_j) == 2:
    #     if 'H' in [at[x] for x in adj_j]:
    #         j_is_oh = True
    if at[i] == 'C':
        # 羟基内部的要算
        # O===C_i---O_j---H
        # if len(adj_i) == 3 and j_is_oh:
        #     for x in adj_i:
        #         if (at[x] == 'O') and (adj_i[x]['type'] == BondType(2)):
        #             return False
        if len(adj_i) == 2:
            k = [x for x in adj_i if x != j][0]
            if adj_i[k]['type'] == 3:
                # C_i ### N
                if at[k] == 'N':
                    return False
                # C_i ### C --- H
                if at[k] == 'C':
                    if 'H' in [at[x] for x in adj[k]]:
                        return False
    return True
