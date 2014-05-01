"""A class representing a chemical substance."""


import numpy
from rdkit import Chem


class Chemical:
    """A class representing a chemical substance."""
    multipliers = {'lynch-willet': 2, 'funatsu': 4}

    def __init__(self, smiles, ec_type='funatsu'):
        """Initializes a chemical.

        Parameters
        ----------
        smiles : string
            A representation of a chemical compound in Daylight SMILES
            notation.

        Raises
        ------
        ValueError
            Exception is thrown if the input SMILES in invalid.
        """
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol is None:
            raise ValueError('Invalid chemical SMILES.')

        # To ensure we are always dealing with a canonical SMILES, we are
        # getting it by converting the molecule back, instead of just using the
        # supplied string.
        self.smiles = Chem.MolToSmiles(self.mol)

        # Save the values of initial indices (within the molecule) on
        # respective atoms; mainly for debugging purposes. Note that those
        # values are deliberately off by one.
        for a in self.mol.GetAtoms():
            a.SetProp('initIdx', str(a.GetIdx() + 1))

        # Get underlying graph matrices.
        self.adjacency_matrix = Chem.GetAdjacencyMatrix(self.mol)
        self.distance_matrix = Chem.GetDistanceMatrix(self.mol)
        self.connectivity_matrix = self.get_acm(self.mol)

        # Attributes for handling extended connectivity (EC) indices.
        # Variable 'ec_order' is rather self-explanatory. It is worth noting
        # though that the element of 'ec_indices' should always correspond
        # to respective atoms, i.e. first element should hold EC index of
        # the first atom; the same should hold for the second and so on.
        self.ec_order = None
        self.ec_indices = []
        self.multi = Chemical.multipliers.get(ec_type, 2)

    def calc_init_ecs(self, index_type='funatsu'):
        """Calculates initial EC indices.

        Parameters
        ----------
        index_type : string
            Initial values of EC indices can be calculated using three
            different methods called 'funatsu', 'shelley', and 'morgan'.
            Defaults to 'funatsu'.

        Returns
        -------
        indices : tuple
            A sequence in which the $i$-th element represent the EC
            index of the corresponding atom.

            .. warning::
            It does *not* update EC values on atoms.
        """
        methods = {
            'funatsu': self.get_funatsu_identifier,
            'morgan': self.get_morgan_identifier,
            'shelley': self.get_shelley_identifier
        }
        get_identifier = methods.get(index_type, self.get_funatsu_identifier)
        return tuple([get_identifier(a) for a in self.mol.GetAtoms()])

    def calc_next_ecs(self):
        """Calculates next EC indices.

        Function calculates next order EC indices for molecule's atom.
        The indices are incremented according to the formula

        .. math::
        EC_{i}^{n} = m * EC_{i}^{n - 1} + \sum_{j} EC_{j}^{n - 1},

        where :math:`m` is method dependent multiplier (4 for 'funatsu',
        2 otherwise) and the summation goes over atoms adjacent to  atom
        $i$.

        Returns
        -------
        indices : tuple
            A sequence in which the $i$-th element represent the EC index
            of the corresponding atom.

            .. warning::
            It does NOT updates their values on atoms.
        """

        # If EC indices were not initialized, do it now.
        if self.ec_order is None:
            self.ec_indices = self.calc_init_ecs()
            self.ec_order = 1

        # Calculate and return higher order EC indices using Lynch-Willett
        # formula.
        ecs = numpy.array(self.ec_indices)
        return tuple(self.multi * ecs + numpy.dot(self.adjacency_matrix, ecs))

    def clear_ecs(self):
        """Clears EC indices."""
        self.ec_order = None
        self.ec_indices = []

    def update_ecs(self, ecs):
        """Updates EC indices.

        Parameters
        ----------
        ecs : sequence
            A sequence of integer representing atoms' EC indices.
        """
        if self.ec_order is None:
            self.ec_order = 1
        else:
            self.ec_order += 1
        self.ec_indices = ecs

    def find_ec_mcs(self, index, radius):
        """Finds all atoms belonging to a given EC-MCS.

        The atom's :math:`n`-order EC index may be treated as a hash of a
        circular substructure with a radius :math:`r \equiv (n - 2)` bonds.

        Parameters
        ----------
        index : integer
            An index of an atom for which EC-MCS will be determined

        radius : integer
            Match radius.

        Returns
        -------
        indices : list
            Indices of all atoms lying *within* match radius, i.e. atoms
            separated by at most :math:`r - 1` bonds from the center atom.
        """
        return [idx for idx, dist in enumerate(self.distance_matrix[index])
                if dist < radius]

    def remove_atoms(self, indices):
        """Removes atoms with given indices."""

        # Every time an atom with is removed from a molecule, atom indices
        # higher than the deleted atom's index are adjusted. Thus, we are
        # removing them in reversed order, i.e. starting from the atom with
        # highest.
        e = Chem.EditableMol(self.mol)
        for i in sorted(indices, reverse=True):
            e.RemoveAtom(i)
        self.mol = e.GetMol()
        self.adjacency_matrix = Chem.GetAdjacencyMatrix(self.mol)
        self.distance_matrix = Chem.GetDistanceMatrix(self.mol)
        self.smiles = Chem.MolToSmiles(self.mol)

    def get_atom_classes(self, atom_indices):
        """Finds out atoms' classes.

        Method finds out the classes of atoms (type and bond pattern)
        from the list of indices supplied as its input argument. Currently,
        it determines an atom's class using method described by Shelley and
        Munk.

        Parameters
        ----------
        atom_indices : sequence
            Indices of atoms that classes need to be determined.

        Returns
        -------
        classes : set
            A set of integers representing atoms' classes.
        """
        return set([self.get_shelley_identifier(self.mol.GetAtomWithIdx(i))
                    for i in atom_indices])

    @staticmethod
    def get_funatsu_identifier(atom):
        """Calculates an initial EC index of a given atom.

        After Funatsu et al. in Tetrahedron Computer Methodology 1,
        53-69 (1988), the values are calculated according to formula

        .. math::
        EC^{1}_{a_{i}} = 10 Z_{a_{i}} + \deg(a_{i})

        where :math:`Z_{a_{i}}` is the atomic number of :math:`i`-th atom.

        Parameters
        ----------
        atom : Atom
            An atom for which calculate index.

        Returns
        -------
        index : integer
            An integer representing initial EC index for a given atom.
        """
        return 10 * atom.GetAtomicNum() + len(atom.GetNeighbors())

    @staticmethod
    def get_morgan_identifier(atom):
        """Calculates an initial EC index of a given atom.

        After Morgan in J. Chem. Doc. 5, 107 (1965), the value is just number
        of nonhydrogen neighbours.

        Parameters
        ----------
        atom : Atom
            An atom for which calculate index.

        Returns
        -------
        index : integer
            An integer representing initial EC index for a given atom.
        """
        return len(atom.GetNeighbors())

    @staticmethod
    def get_shelley_identifier(atom):
        """Calculates and returns an initial EC index of a given atom.

        After Shelley and Munk in J. Chem. Inf. Com. Sci 17, 110 (1977),
        the values are two-digit integer, the most significant of which
        specifies number of covalent (two-electrons) bonds by which it joins
        nonhydrogen atoms, and the least significant of which designates
        atom elemental type (i.e. :math:`C = 2`, :math:`N = 3`, :math:`O = 4`,
        etc.)

        Parameters
        ----------
        atom : Atom
            An atom for which calculate index.

        Returns
        -------
        index : integer
            An integer representing initial EC index for a given atom.
        """

        # Get atom's type. If not defined, revert to atomic number.
        organic_subset = {
            'B': 1, 'C': 2, 'N': 3, 'O': 4, 'P': 5,
            'S': 6, 'F': 7, 'Cl': 8, 'Br': 9, 'I': 10
        }
        idx = organic_subset.get(atom.GetSymbol(), atom.GetAtomicNum())
        valence = sum([int(b.GetValenceContrib(atom))
                       for b in atom.GetBonds()])
        return 10 * valence + idx

    @staticmethod
    def get_acm(mol):
        """Calculates atom connectivity matrix of a molecule.

        Atom connectivity matrix is a weighted adjacency matrix which
        diagonal elements indicates atom's atomic number and non-zero
        off-diagonal elements denotes bond order.

        Parameters
        ----------
        mol : Molecule
            A molecule for which

        Returns
        -------
        matrix : numpy array
            An array representing atom connectivity matrix.
        """
        matrix = numpy.zeros((len(mol.GetAtoms()), len(mol.GetAtoms())))
        for i, a in enumerate(mol.GetAtoms()):
            matrix[i][i] = a.GetAtomicNum()
            for j in [a.GetIdx() for a in a.GetNeighbors()]:
                matrix[i][j] = mol.GetBondBetweenAtoms(i, j).GetBondType()
        return matrix


if __name__ == '__main__':
    smi = 'CC1CCCC2CCCC(C)C12'
    chem = Chemical(smi)
    print chem.smiles
    print 'First order ECs:', chem.calc_init_ecs(index_type='morgan')
    print 'Second order ECs:', chem.calc_next_ecs(index_type='morgan')
