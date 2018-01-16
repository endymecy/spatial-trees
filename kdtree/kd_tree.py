# -*- coding: utf-8 -*-

from collections import deque
from kdtree.kd_node import KDNode


def create(point_list=None, dimensions=None, axis=0, sel_axis=None):
    """ Creates a kd-tree from a list of points
        All points in the list must be of the same dimensionality.
        If no point_list is given, an empty tree is created. The number of
        dimensions has to be given instead.
        If both a point_list and dimensions are given, the numbers must agree.
        Axis is the axis on which the root-node should split.
        sel_axis(axis) is used when creating sub_nodes of a node. It receives the
        axis of the parent node and returns the axis of the child node. """
    if not point_list and not dimensions:
        raise ValueError('either point_list or dimensions must be provided')
    elif point_list:
        dimensions = check_dimensionality(point_list, dimensions)

    # 默认是顺序选择，也可以使用方差来选择
    sel_axis = sel_axis or (lambda prev_axis: (prev_axis + 1) % dimensions)

    # 返回空的树
    if not point_list:
        return KDNode(sel_axis=sel_axis, axis=axis, dimensions=dimensions)

    # 对list进行排序并选择中位数作为切分点
    point_list = list(point_list)
    point_list.sort(key=lambda point: point[axis])
    median = len(point_list) // 2

    loc = point_list[median]
    left = create(point_list[:median], dimensions, sel_axis(axis))
    right = create(point_list[median + 1:], dimensions, sel_axis(axis))
    return KDNode(loc, left, right, axis=axis, sel_axis=sel_axis, dimensions=dimensions)


def check_dimensionality(point_list, dimensions=None):
    dimensions = dimensions or len(point_list[0])
    for p in point_list:
        if len(p) != dimensions:
            raise ValueError('All Points in the point_list must have the same dimensionality')

    return dimensions


def level_order(tree, include_all=False):
    """ Returns an iterator over the tree in level-order
    If include_all is set to True, empty parts of the tree are filled
    with dummy entries and the iterator becomes infinite. """

    q = deque()
    q.append(tree)
    while q:
        node = q.popleft()
        yield node

        if include_all or node.left:
            q.append(node.left or node.__class__())

        if include_all or node.right:
            q.append(node.right or node.__class__())
