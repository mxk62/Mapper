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
        self.d = Chem.GetDistanceMatrix(self.mol)

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
        for a in self.mol.GetAtoms():
            ecn = 2 * int(a.GetProp('EC')) + \
                sum(int(n.GetProp('EC')) for n in a.GetNeighbors())
            a.SetProp('EC', str(ecn))
        self.ec_order += 1
        return tuple(int(a.GetProp('EC')) for a in self.mol.GetAtoms())

    def init_ecs(self):
        """Sets and returns a tuple representing initial EC indices.

        Function initializes EC indices on the molecule. Their are stored as
        the 'EC' atom property. Right now, the initial values are just the
        number of non-hydrogen bonds attached to a given atom.

        It returns a tuple in which the $i$-th element represent the EC
        index of the corresponding atom.
        """
        for a in self.mol.GetAtoms():
            ec = len(a.GetNeighbors())
            a.SetProp('EC', str(ec))
        self.ec_order = 0
        return tuple(int(a.GetProp('EC')) for a in self.mol.GetAtoms())

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

    def find_ec_mcs(self, index):
        """Returns indices of all atoms belonging to a given EC-MCS.

        The atom's $n$-order EC index may be treated as a hash of a circular
        substructure with a radius $(n - 1)$ bonds. Function returns indices of
        ALL atoms lying within.
        """
        return [idx for idx, dist in enumerate(self.d[index])
                if dist < self.ec_order - 1]


if __name__ == '__main__':
    smi = 'CC1CCCC2CCCC(C)C12'
    chem = Chemical(smi)
    chem.init_ecs()
    for a in chem.mol.GetAtoms():
        print a.GetSymbol(), a.GetIdx(), a.GetProp('EC')
    chem.increment_ecs()
    for a in chem.mol.GetAtoms():
        print a.GetSymbol(), a.GetIdx(), a.GetProp('EC')
