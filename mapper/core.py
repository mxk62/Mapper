"""Class representing reaction core."""


import itertools
from rdkit import Chem
from rdkit.Chem import SanitizeFlags
import mapper as mp


class Core:
    """Represents a reaction core."""

    def __init__(self, smarts):
        self.smarts = smarts

        frags = [frag.split('.') for frag in smarts.split('>>')]

        self.retrons, self.synthons = \
            [[Chem.MolFromSmarts(smi) for smi in smis] for smis in frags]

        self.reactants, self.products = \
            [[Chem.MolFromSmiles(smi, sanitize=False) for smi in smis]
             for smis in frags]
        opts = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        for m in self.reactants + self.products:
            Chem.SanitizeMol(m, sanitizeOps=opts)

    def __eq__(self, other):
        return self.smiles == self.other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.smarts)

    def does_contain(self, other):
        """Returns True if contains core.

        Parameters
        ----------
        other : another reaction core

        Examples
        --------
        Core '[C:1][O:2]>>[C:1]=[O:2]' is considered to contain core
        '[C:2]>>[C:2]Br' as the latter core's retron, '[C:2]', is
        a substructure of '[C:1][O:2]'.

        >>> core = mp.Core('[C:1][O:2]>>[C:1]=[O:2]')
        >>> other = mp.Core('[C:2]>>[C:2]Br')
        >>> core.does_contain(other)
        True

        Conversely, it does not contain the core '[S:2]>>O=[S:2]=O' as the
        the other core's retron cannot be a substructure of '[C:1][O:2]'.

        >>> core = mp.Core('[C:1][O:2]>>[C:1]=[O:2]')
        >>> other = mp.Core('[S:2]>>O=[S:2]=O')
        >>> core.does_contain(other)
        False

        Core having multiple retrons is considered to contain other
        multiple-retron core if and only if:
          1. number of retrons agrees between the cores;
          2. for *any* ordering of the latter core's retrons, *simultanously*
             all of them are substructures of the former core's retrons.
        Thus

        >>> core = mp.Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
        >>> other = mp.Core('[N:3].[C:2]O>>[C:2][N:3]')
        >>> core.does_contain(other)
        True

        but

        >>> core = mp.Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
        >>> other = mp.Core('[C:1]=O>>[C:1]=C')
        >>> core.does_contain(other)
        False
        """
        if len(self.reactants) != len(other.retrons):
            return False
        mols = itertools.permutations(self.reactants)
        patterns = itertools.permutations(other.retrons)
        for mols, patterns in zip(mols, patterns):
            if all(m.HasSubstruct(p) for m, p in zip(mols, patterns)):
                return True
        return False