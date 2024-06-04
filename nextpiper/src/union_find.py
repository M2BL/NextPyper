# =======================================================================================
#               CLASS
# =======================================================================================


class UnionFind:
    """Union-find data structure with path compression and Union-by-rank.
    getComponents return a list of sets, each set instance represents
    a component.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}
        self.components = {}

    def __getitem__(self, object):
        """
        Find and return the name of the set containing the object.
        """
        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            self.components[object] = [object]
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]
        for ancestor in path:
            self.parents[ancestor] = root
            if ancestor != root:
                self.components[root].append(ancestor)
                try:
                    del self.components[ancestor]
                except:
                    continue
        return root

    def union(self, objects):
        """
        Find the sets containing the objects and merge them all.
        """
        roots = [self[x] for x in objects]  # uses__getitem__
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest
                self.components[heaviest].extend(self.components[r])
                del self.components[r]

    def get_components(self):
        return [set(x) for x in self.components.values()]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
