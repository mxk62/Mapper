import numpy
from rdkit import Chem


class Chemical:
    """A class representing a chemical substance."""

    def __init__(self, smiles):
        try:
            self.mol = Chem.MolFromSmiles(smiles)
        except:
            raise ValueError('Invalid chemical SMILES.')
        self.smiles = Chem.MolToSmiles(self.mol)
        self.ec_order = None
        self.adjacency_matrix = Chem.GetAdjacencyMatrix(self.mol)
        self.distance_matrix = Chem.GetDistanceMatrix(self.mol)

    def increment_ecs(self):
        """Sets and returns a tuple representing atom EC indices.

        Function calculates next order EC index for each atom in the molecule,
        replacing the old values. The indices are increment using Lynch-Willet
        formula, i.e.
        \[
            EC_{i}^{n} = 2 * EC_{i}^{n - 1} + \sum_{j} EC_{j}^{n - 1},
        \]
        where the summation goes over atoms adjacent to atom $i$.

        It returns a tuple in which the $i$-th element represent the EC index
        of the corresponding atom.
        """

        # If EC indices were not initialized, do it now.
        if self.ec_order is None:
            self.init_ecs()

        # Get current EC indices.
        ecs = numpy.array([int(a.GetProp('EC')) for a in self.mol.GetAtoms()])

        # Calculate higher order EC indices using Lynch-Willett formula.
        ecs = 2 * ecs + numpy.dot(self.adjacency_matrix, ecs)

        # Update atom properties and order value.
        for idx, ecn in enumerate(ecs):
            self.mol.GetAtomWithIdx(idx).SetProp('EC', str(ecn))
        self.ec_order += 1
        return tuple(ecs)

    def init_ecs(self):
        """Sets and returns a tuple representing initial EC indices.

        Function initializes EC indices on the molecule. After Funatsu et al.
        in Tetrahedron Computer Methodology 1, 53-69 (1988), the initial values
        are calculated according to formula
        \[
            EC^{1}_{a_{i}} = 10 Z_{a_{i}} + \deg(a_{i})
        \]
        where $Z_{a_{i}}$ is the atomic number of $i$-th atom.

        Their are stored as the 'EC' atom property. Function returns a
        tuple in which the $i$-th element represent the EC index of the
        corresponding atom.
        """
        for a in self.mol.GetAtoms():
            a.SetProp('EC', str(10 * a.GetAtomicNum() + len(a.GetNeighbors())))
        self.ec_order = 0
        return tuple(int(a.GetProp('EC')) for a in self.mol.GetAtoms())

    def clear_ecs(self):
        """Clear EC indices from molecule's atoms.

        Function literally removes EC indices from molecule's atoms, i.e.
        the property is removed along with the value it holds.
        """
        for a in self.mol.GetAtoms():
            a.ClearProp('EC')
        self.ec_order = None

    def remove_atoms(self, indices):
        """Removes atoms with given indices.

        Every time an atom with is removed from a molecule, atom indices
        higher than the deleted atom's index are adjusted. Thus, we are
        removing them in reversed order, i.e. starting from the atom with
        highest.
        """
        e = Chem.EditableMol(self.mol)
        for i in sorted(indices, reverse=True):
            e.RemoveAtom(i)
        self.mol = e.GetMol()
        self.adjacency_matrix = Chem.GetAdjacencyMatrix(self.mol)
        self.distance_matrix = Chem.GetDistanceMatrix(self.mol)
        self.smiles = Chem.MolToSmiles(self.mol)

    def find_ec_mcs(self, index):
        """Returns indices of all atoms belonging to a given EC-MCS.

        The atom's $n$-order EC index may be treated as a hash of a circular
        substructure with a radius $(n - 1)$ bonds. Function returns indices of
        ALL atoms lying within.
        """
        return [idx for idx, dist in enumerate(self.distance_matrix[index])
                if dist < self.ec_order]


if __name__ == '__main__':
    smi = 'CC1CCCC2CCCC(C)C12'
    chem = Chemical(smi)
    chem.init_ecs()
    for a in chem.mol.GetAtoms():
        print a.GetSymbol(), a.GetIdx(), a.GetProp('EC')
    chem.increment_ecs()
    for a in chem.mol.GetAtoms():
        print a.GetSymbol(), a.GetIdx(), a.GetProp('EC')
