# -*- coding: utf-8 -*-

from . import kd_tree
import math
import itertools
import heapq


class Node(object):
    def __init__(self, data=None, left=None, right=None):
        self.data = data
        self.left = left
        self.right = right

    @property
    def is_leaf(self):
        return (not self.data) or (all(not bool(c) for c, _ in self.children))

    def pre_order(self):
        if not self:
            return

        yield self

        if self.left:
            for x in self.left.pre_order():
                yield x

        if self.right:
            for x in self.right.pre_order():
                yield x

    def in_order(self):
        if not self:
            return

        if self.left:
            for x in self.left.in_order():
                yield x

        yield self

        if self.right:
            for x in self.right.in_order():
                yield x

    def post_order(self):
        if not self:
            return

        if self.left:
            for x in self.left.post_order():
                yield x

        if self.right:
            for x in self.right.post_order():
                yield x

        yield self

    @property
    def children(self):
        """
        Returns an iterator for the non-empty children of the Node
        The children are returned as (Node, pos) tuples where pos is 0 for the
        left sub_node and 1 for the right_node.
        """
        if self.left and self.left.data is not None:
            yield self.left, 0
        if self.right and self.right.data is not None:
            yield self.right, 1

    def set_child(self, index, child):
        if index == 0:
            self.left = child
        else:
            self.right = child

    def height(self):
        min_height = int(bool(self))
        return max([min_height] + [c.height() + 1 for c, p in self.children])

    def get_child_pos(self, child):
        for c, pos in self.children:
            if child == c:
                return pos

    def __eq__(self, other):
        if isinstance(other, tuple):
            return self.data == other
        else:
            return self.data == other.data

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return '<%(cls)s - %(data)s>' % \
               dict(cls=self.__class__.__name__, data=repr(self.data))

    def __nonzero__(self):
        return self.data is not None

    __bool__ = __nonzero__


class KDNode(Node):
    """ A Node that contains kd-tree specific data and methods """

    def __init__(self, data=None, left=None, right=None, axis=None,
                 sel_axis=None, dimensions=None):
        """ Creates a new node for a kd-tree
            If the node will be used within a tree, the axis and the sel_axis
            function should be supplied.
            sel_axis(axis) is used when creating subnodes of the current node. It
            receives the axis of the parent node and returns the axis of the child
            node. """
        super(KDNode, self).__init__(data, left, right)
        self.axis = axis  # 选中的特征列
        self.sel_axis = sel_axis  # 如何选择特征列
        self.dimensions = dimensions

    def add(self, point):
        current = self
        while True:
            kd_tree.check_dimensionality([point], dimensions=current.dimensions)

            # 空的子节点，直接添加
            if current.data is None:
                current.data = point
                return current

            # split on self.axis, recurse either left or right
            # 判断添加点的current.axis数据是否小于当前节点的current.axis数据
            if point[current.axis] < current.data[current.axis]:
                if current.left is None:
                    current.left = current.create_subnode(point)
                    return current.left
                else:
                    current = current.left
            else:
                if current.right is None:
                    current.right = current.create_subnode(point)
                    return current.right
                else:
                    current = current.right

    def create_subnode(self, data):
        """ Creates a subnode for the current node """
        return self.__class__(data,
                              axis=self.sel_axis(self.axis),
                              sel_axis=self.sel_axis,
                              dimensions=self.dimensions)

    def find_replacement(self):
        if self.right:
            child, parent = self.right.extreme_child(min, self.axis)
        else:
            child, parent = self.left.extreme_child(max, self.axis)

        return child, parent if parent is not None else self

    def should_remove(self, point, node):
        """ checks if self's point (and maybe identity) matches """
        if not self.data == point:
            return False

        return (node is None) or (node is self)

    def remove(self, point, node=None):
        if not self:
            return

        # Recursion has reached the node to be deleted
        if self.should_remove(point, node):
            return self._remove(point)

        # Remove direct subnode
        if self.left and self.left.should_remove(point, node):
            self.left = self.left._remove(point)
        elif self.right and self.right.should_remove(point, node):
            self.right = self.right._remove(point)

        # Recurse to subtrees
        if point[self.axis] <= self.data[self.axis]:
            if self.left:
                self.left = self.left.remove(point, node)

        if point[self.axis] >= self.data[self.axis]:
            if self.right:
                self.right = self.right.remove(point, node)

        return self

    def _remove(self, point):

        if self.is_leaf:
            self.data = None
            return self

        root, max_p = self.find_replacement()

        # self and root swap positions
        tmp_l, tmp_r = self.left, self.right
        self.left, self.right = root.left, root.right
        root.left, root.right = tmp_l if tmp_l is not root else self, tmp_r if tmp_r is not root else self
        self.axis, root.axis = root.axis, self.axis

        # Special-case if we have not chosen a direct child as the replacement
        if max_p is not self:
            pos = max_p.get_child_pos(root)
            max_p.set_child(pos, self)
            max_p.remove(point, self)

        else:
            root.remove(point, self)

        return root

    @property
    def is_balanced(self):
        """ Returns True if the (sub)tree is balanced
        The tree is balanced if the heights of both subtrees differ at most by
        1 """

        left_height = self.left.height() if self.left else 0
        right_height = self.right.height() if self.right else 0

        if abs(left_height - right_height) > 1:
            return False

        return all(c.is_balanced for c, _ in self.children)

    def rebalance(self):
        """
        Returns the (possibly new) root of the rebalanced tree
        """

        return kd_tree.create([x.data for x in self.in_order()])

    def axis_dist(self, point, axis):
        """
        Squared distance at the given axis between
        the current Node and the given point
        """
        return math.pow(self.data[axis] - point[axis], 2)

    def dist(self, point):
        """
        Squared distance between the current Node
        and the given point
        """
        r = range(self.dimensions)
        return sum([self.axis_dist(point, i) for i in r])

    def search_knn(self, point, k, dist=None):
        """ Return the k nearest neighbors of point and their distances
            point must be an actual point, not a node.
            k is the number of results to return. The actual results can be less
            (if there aren't more nodes to return) or more in case of equal
            distances.
            dist is a distance function, expecting two points and returning a
            distance value. Distance values can be any comparable type.
            The result is an ordered list of (node, distance) tuples.
        """
        if k < 1:
            raise ValueError("k must be greater than 0.")

        if dist is None:
            get_dist = lambda n: n.dist(point)
        else:
            get_dist = lambda n: dist(n.data, point)

        results = []

        self._search_node(point, k, results, get_dist, itertools.count())

        return [(node, -d) for d, _, node in sorted(results, reverse=True)]

    def _search_node(self, point, k, results, get_dist, counter):
        if not self:
            return

        node_dist = get_dist(self)  # 计算当前点和point的距离

        # Add current node to the priority queue if it closer than
        # at least one point in the queue.
        #
        # If the heap is at its capacity, we need to check if the
        # current node is closer than the current farthest node, and if
        # so, replace it.

        item = (-node_dist, next(counter), self)
        if len(results) >= k:
            if -node_dist > results[0][0]:
                heapq.heapreplace(results, item)
        else:
            heapq.heappush(results, item)
        # get the splitting plane
        split_plane = self.data[self.axis]  # 中位数即切分点
        # get the squared distance between the point and the splitting plane
        # (squared since all distances are squared).
        plane_dist = point[self.axis] - split_plane
        plane_dist2 = plane_dist * plane_dist

        # Search the side of the splitting plane that the point is in
        if point[self.axis] < split_plane:
            if self.left is not None:
                self.left._search_node(point, k, results, get_dist, counter)
        else:
            if self.right is not None:
                self.right._search_node(point, k, results, get_dist, counter)

        # Search the other side of the splitting plane if it may contain
        # points closer than the farthest point in the current results.
        if -plane_dist2 > results[0][0] or len(results) < k:  # 回溯查找
            if point[self.axis] < self.data[self.axis]:
                if self.right is not None:
                    self.right._search_node(point, k, results, get_dist, counter)
            else:
                if self.left is not None:
                    self.left._search_node(point, k, results, get_dist, counter)

    def search_nn(self, point, dist=None):
        """
        Search the nearest node of the given point
        point must be an actual point, not a node. The nearest node to the
        point is returned. If a location of an actual node is used, the Node
        with this location will be returned (not its neighbor).
        dist is a distance function, expecting two points and returning a
        distance value. Distance values can be any comparable type.
        The result is a (node, distance) tuple.
        """

        return next(iter(self.search_knn(point, 1, dist)), None)

    def _search_nn_dist(self, point, dist, results, get_dist):
        if not self:
            return

        node_dist = get_dist(self)

        if node_dist < dist:
            results.append(self.data)

        # get the splitting plane
        split_plane = self.data[self.axis]

        # Search the side of the splitting plane that the point is in
        if point[self.axis] <= split_plane + dist:
            if self.left is not None:
                self.left._search_nn_dist(point, dist, results, get_dist)
        if point[self.axis] >= split_plane - dist:
            if self.right is not None:
                self.right._search_nn_dist(point, dist, results, get_dist)

    def search_nn_dist(self, point, distance, best=None):
        """
        Search the n nearest nodes of the given point which are within given
        distance
        point must be a location, not a node. A list containing the n nearest
        nodes to the point within the distance will be returned.
        """

        results = []
        get_dist = lambda n: n.dist(point)

        self._search_nn_dist(point, distance, results, get_dist)
        return results

    def is_valid(self):
        """ Checks recursively if the tree is valid
        It is valid if each node splits correctly """

        if not self:
            return True

        if self.left and self.data[self.axis] < self.left.data[self.axis]:
            return False

        if self.right and self.data[self.axis] > self.right.data[self.axis]:
            return False

        return all(c.is_valid() for c, _ in self.children) or self.is_leaf

    def extreme_child(self, sel_func, axis):
        """ Returns a child of the subtree and its parent
        The child is selected by sel_func which is either min or max
        (or a different function with similar semantics). """

        max_key = lambda child_parent: child_parent[0].data[axis]

        # we don't know our parent, so we include None
        me = [(self, None)] if self else []

        child_max = [c.extreme_child(sel_func, axis) for c, _ in self.children]
        # insert self for unknown parents
        child_max = [(c, p if p is not None else self) for c, p in child_max]

        candidates = me + child_max

        if not candidates:
            return None, None

        return sel_func(candidates, key=max_key)




