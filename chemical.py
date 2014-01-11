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
        self.adjacency_matrix = Chem.GetAdjacencyMatrix(self.mol)
        self.distance_matrix = Chem.GetDistanceMatrix(self.mol)

        # Attributes for handling extended connectivity (EC) indices.
        # Variable 'ec_order' is rather self-explanatory. It is worth noting
        # though that the element of 'ec_indices' should always correspond
        # to respective atoms, i.e. first element should hold EC index of
        # the first atom; the same should hold for the second and so on.
        self.ec_order = None
        self.ec_indices = []

    def calc_init_ecs(self, index_type='funatsu'):
        """Calculates and returns a tuple representing initial EC indices.

        Function returns a tuple in which the $i$-th element represent the EC
        index of the corresponding atom. It does NOT update their values on
        atoms.
        """
        methods = {
            'funatsu': self.get_funatsu_identifier,
            'morgan': self.get_morgan_identifier,
            'shelley': self.get_shelley_identifier
        }
        get_identifier = methods.get(index_type, self.get_funatsu_identifier)
        return tuple([get_identifier(a) for a in self.mol.GetAtoms()])

    @staticmethod
    def get_funatsu_identifier(atom):
        """Calculates and returns an initial EC index of a given atom.

        After Funatsu et al. in Tetrahedron Computer Methodology 1,
        53-69 (1988), the values are calculated according to formula
        \[
            EC^{1}_{a_{i}} = 10 Z_{a_{i}} + \deg(a_{i})
        \]
        where $Z_{a_{i}}$ is the atomic number of $i$-th atom.
        """
        return 10 * atom.GetAtomicNum() + len(atom.GetNeighbors())

    @staticmethod
    def get_morgan_identifier(atom):
        """Calculates and returns an initial EC index of a given atom.

        After Morgan ... the values is just number of nonhydrogen neighbours.
        """
        return len(atom.GetNeighbors())

    @staticmethod
    def get_shelley_identifier(atom):
        """Calculates and returns an initial EC index of a given atom.

        After Shelley and Munk in J. Chem. Inf. Com. Sci 17, 110 (1977),
        the values are calculated as a sum of atom's elemental type
        (C = 2, N = 3, O = 4, etc.) times 10 and number of covalent
        (two-electrons) bonds by which it joins nonhydrogen atoms.
        """

        # Get atom's type. If not defined, revert to atomic number.
        organic_subset = {
            'B': 1, 'C': 2, 'N': 3, 'O': 4, 'P': 5,
            'S': 6, 'F': 7, 'Cl': 8, 'Br': 9, 'I': 10
        }
        idx = organic_subset.get(atom.GetSymbol(), atom.GetAtomicNum())
        valence = sum([b.GetValenceContrib(atom) for b in atom.GetBonds()])
        return 10 * idx + int(valence)

    def calc_next_ecs(self):
        """Calculates and returns a tuple representing next EC indices.

        Function calculates next order EC indices for molecule's atom.
        The indices are incremented according to Lynch-Willet formula, i.e.
        \[
            EC_{i}^{n} = 2 * EC_{i}^{n - 1} + \sum_{j} EC_{j}^{n - 1},
        \]
        where the summation goes over atoms adjacent to atom $i$.

        It returns a tuple in which the $i$-th element represent the EC index
        of the corresponding atom. It does NOT updates their values on atoms.
        """

        # If EC indices were not initialized, do it now.
        if self.ec_order is None:
            self.ec_indices = self.calc_init_ecs()
            self.ec_order = 1

        # Calculate and return higher order EC indices using Lynch-Willett
        # formula.
        ecs = numpy.array(self.ec_indices)
        return tuple(2 * ecs + numpy.dot(self.adjacency_matrix, ecs))

    def clear_ecs(self):
        """Clear EC indices from molecule's atoms."""
        self.ec_order = None
        self.ec_indices = []

    def update_ecs(self, ecs):
        """Updates EC indices on atoms and adjust their order."""
        if self.ec_order is None:
            self.ec_order = 1
        else:
            self.ec_order += 1
        self.ec_indices = ecs

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

    def find_ec_mcs(self, index, radius):
        """Returns indices of all atoms belonging to a given EC-MCS.

        The atom's $n$-order EC index may be treated as a hash of a circular
        substructure with a radius $(n - 2)$ bonds. Function returns indices of
        ALL atoms lying within.
        """
        return [idx for idx, dist in enumerate(self.distance_matrix[index])
                if dist <= radius]


if __name__ == '__main__':
    smi = 'CC1CCCC2CCCC(C)C12'
    chem = Chemical(smi)
    chem.calc_init_ecs()
    chem.calc_next_ecs()
