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


if __name__ == 'main':
    smi = 'C'
    chem = Chemical(smi)
    for a in chem.mol.GetAtoms():
        print a.GetSymbol(), a.GetIdx(), a.GetProp('EC')
