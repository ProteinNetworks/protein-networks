"""
In order to compare two proteins or protein fragments, we should have some
network-based comparison measure.

Use the hierarchical community structure to compare the networks.


Steps:
- Convert the dendrogram into a tree
    - How? Weighting?
- Get some distance measure
    - Pre-existing weighted tree distances or roll-your-own, idk.
???
- Check TM-score/other measures correlation

Need generated-tree viewer as I go.
"""

from .partition import Partition


class Tree:
    """Holds the various tree things."""
    def __init__(self, inputPartition: Partition):
        """Generate a tree from a given partition."""
        pass

    