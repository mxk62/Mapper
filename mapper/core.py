"""Class representing reaction core."""


import itertools
from rdkit import Chem
from rdkit.Chem import SanitizeFlags
import mapper as mp


class Core:
    """Represents a reaction core.

    Reaction center indicates bond(s) broken or formed during the reaction.
    Reaction core represents substructure (bonds *and* atoms) associated with
    the reaction center.
    """

    def __init__(self, smarts):
        """Initializes a reaction core.

        Parameters
        ----------
        smarts : reaction core SMARTS
            Reaction core encoded using Daylight SMARTS notation.

        Examples
        --------
        >>> c = mp.Core('[C:1][O:2]>>[C:1]=[O:2]')

        Cores represent *patterns* but can be also view as *objects*
        (molecular fragments), thus

        >>> all(m.HasSubstruct(p) for m, p in zip(c.reactants, c.retrons))
        True
        """
        self.smarts = smarts

        frags = [frag.split('.') for frag in smarts.split('>>')]

        self.retrons, self.synthons = \
            [[Chem.MolFromSmarts(smi) for smi in smis] for smis in frags]
        if None in self.retrons + self.synthons:
            raise ValueError('erroneous core SMARTS')

        reactants, products = \
            [[Chem.MolFromSmiles(smi, sanitize=False) for smi in smis]
             for smis in frags]
        self.reactants = self._strip(reactants)
        self.products = self._strip(products)
        opts = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        for m in self.reactants + self.products:
            Chem.SanitizeMol(m, sanitizeOps=opts)
        self.smiles = '>>'.join([
            '.'.join([Chem.MolToSmiles(m) for m in self.reactants]),
            '.'.join([Chem.MolToSmiles(m) for m in self.products])
        ])

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.smiles)

    def does_contain(self, other):
        """Returns True if

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
        for mols in itertools.permutations(self.reactants):
            if all(m.HasSubstructMatch(p)
                   for m, p in zip(mols, other.retrons)):
                return True
        return False

    def _strip(self, mols):
        """Strips reactants off unwanted parts."""
        return self._strip_map(self._strip_env(mols))

    @staticmethod
    def _strip_env(mols):
        """Removes unmapped atoms."""
        for m in mols:
            indices = [a.GetIdx() for a in m.GetAtoms()
                       if not a.HasProp('molAtomMapNumber')]

            e = Chem.EditableMol(m)
            for i in sorted(indices, reverse=True):
                e.RemoveAtom(i)

            m = e.GetMol()
        return mols

    @staticmethod
    def _strip_map(mols):
        """Removes atoms map numbers."""
        for m in mols:
            [a.ClearProp('molAtomMapNumber')
             for a in m.GetAtoms() if a.HasProp('molAtomMapNumber')]
        return mols
