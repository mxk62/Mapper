"""Class representing a reaction core.

Reaction center indicates bond(s) broken or formed during the reaction.
Reaction core represents substructure (bonds *and* atoms) associated with
the reaction center.
"""

__all__ = ["Core"]


import itertools

from rdkit import Chem
from rdkit.Chem import MCS, SanitizeFlags


class Core:
    """Class representing a reaction core.

    Parameters
    ----------
    smarts : str
        Reaction core encoded using Daylight SMARTS notation.

    Raises
    ------
    ValueError
        Raised if the input SMARTS is invalid.

    Examples
    --------
    >>> c = Core('[C:1][O:2]>>[C:1]=[O:2]')

    Core parts represent *patterns* but can be also view as
    *objects* (molecular fragments), thus

    >>> all(m.HasSubstructMatch(p) for m, p in zip(c.reactants, c.retrons))
    True

    The same holds for synthons too, e.g.

    >>> all(m.HasSubstructMatch(p) for m, p in zip(c.products, c.synthons))
    True

    """

    def __init__(self, smarts: str):
        self.smarts = smarts

        frags = [frag.split('.') for frag in smarts.split('>>')]

        self.retrons, self.synthons = \
            [[Chem.MolFromSmarts(smi) for smi in smis] for smis in frags]
        if None in self.retrons + self.synthons:
            raise ValueError('erroneous core SMARTS')

        reactants, products = \
            [[Chem.MolFromSmiles(smi, sanitize=False) for smi in smis]
             for smis in frags]
        opts = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        for m in reactants + products:
            Chem.SanitizeMol(m, sanitizeOps=opts)
        self.reactants = self._strip_map(reactants)
        self.products = self._strip_map(products)
        self.smiles = '>>'.join([
            '.'.join([Chem.MolToSmiles(m) for m in self.reactants]),
            '.'.join([Chem.MolToSmiles(m) for m in self.products])
        ])

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Core):
            return NotImplemented
        return self.smiles == other.smiles

    def __ne__(self, other: object) -> bool:
        if not isinstance(other, Core):
            return NotImplemented
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash(self.smiles)

    def __str__(self) -> str:
        return self.smarts

    def does_contain(self, other: "Core") -> bool:
        """Returns True if a core's retrons contains retrons of the other one.

        A retron is considered to contain another one if the latter is the
        substructure of the first one.

        Parameters
        ----------
        other : Core
            Another reaction core

        Examples
        --------
        Core ``[C:1][O:2]>>[C:1]=[O:2]`` contain core ``[C:2]>>[C:2]Br`` as
        the latter core's retron, ``[C:2]``, is a substructure of
        ``[C:1][O:2]``.

        >>> core = Core('[C:1][O:2]>>[C:1]=[O:2]')
        >>> other = Core('[C:2]>>[C:2]Br')
        >>> core.does_contain(other)
        True

        Conversely, it does not contain the core ``[S:2]>>O=[S:2]=O`` as the
        retron ``[S:2]`` is not a substructure of ``[C:1][O:2]``.

        >>> core = Core('[C:1][O:2]>>[C:1]=[O:2]')
        >>> other = Core('[S:2]>>O=[S:2]=O')
        >>> core.does_contain(other)
        False

        Core having multiple retrons is considered to contain other
        multiple-retron core if and only if:

        1. number of retrons agrees between the cores;
        2. for *any* ordering of the latter core's retrons, *simultanously*
           all of them are substructures of the former core's retrons.

        Thus

        >>> core = Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
        >>> other = Core('[N:3].[C:2]O>>[C:2][N:3]')
        >>> core.does_contain(other)
        True

        but

        >>> core = Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
        >>> other = Core('[C:1]=O>>[C:1]=C')
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

    def find_distance(self, other: "Core"):
        r"""Finds similarity distance between two molecular fragments.

        Parameters
        ----------
        other : Core
            The other molecule

        Returns
        -------
        distance : float
            Similarity distance between this and other molecule.

        Note
        ----
        Similarity distance of two molecules is a *metric* defined as

        .. math::
           d(m_{1}, m_{2}) = 1 -
             \frac{|\text{MCS}\{m_{1}, m_{2}\}|}{\max(|m_{1}|, |m_{2}|)}

        where :math:`|\ldots|` denotes number of atoms in a given molecule or
        fragment.

        For any molecules :math:`m_{1}`, :math:`m_{2}` and :math:`m_{3}`,
        the following properties hold true:

        1. :math:`0 \leq d(m_{1}, m_{2}) \le 1`,
        2. :math:`d(m_{1}, m_{2}) = 0` if and only if :math:`m_{1}` and
           :math:`m_{2}` are identical,
        3. :math:`d(m_{1}, m_{2}) = d(m_{2}, m_{1})`,
        4. :math:`d(m_{1}, m_{3}) \leq d(m_{1}, m_{2}) + d(m_{2}, m_{3})`
        """
        if len(self.reactants) != len(other.reactants):
            return 1.0

        distances = []
        for mols in itertools.permutations(self.reactants):
            d = 0.0
            for idx, (moli, molj) in enumerate(zip(mols, other.reactants)):
                moli_size = len(moli.GetAtoms())
                molj_size = len(molj.GetAtoms())

                mcs_size = 0
                if moli_size == 1 or molj_size == 1:
                    if (moli.HasSubstructMatch(other.retrons[idx]) or
                            molj.HasSubstructMatch(self.retrons[idx])):
                        mcs_size = 1
                else:
                    mcs = MCS.FindMCS([moli, molj], ringMatchesRingOnly=True)
                    if mcs.numAtoms != -1:
                        mcs_size = mcs.numAtoms

                d += 1.0 - float(mcs_size) / max(moli_size, molj_size)
            distances.append(d / len(mols))
        return min(distances)

    def _strip(self, mols: list[Chem]) -> list[Chem]:
        """Strips reactants off unwanted elements.

        Parameters
        ---------
        mols : list[Chem]
            List of reactants to clean up.

        Returns
        -------
        mols : list[Chem]
            List of reactants with unwanted elements removed.
        """

        # The order of operations is important here as stripping atom map
        # numbers first will lead to total removal of the molecules.
        return self._strip_map(self._strip_env(mols))

    @staticmethod
    def _strip_env(mols: list[Chem]) -> list[Chem]:
        """Removes unmapped atoms.

        Parameters
        ---------
        mols : list[Chem]
            List of molecules to clean up.

        Returns
        -------
        mols : list[Chem]
            List of molecules with unmapped atoms removed.
        """
        for m in mols:
            indices = [a.GetIdx() for a in m.GetAtoms()
                       if not a.HasProp('molAtomMapNumber')]

            e = Chem.EditableMol(m)
            for i in sorted(indices, reverse=True):
                e.RemoveAtom(i)

            m = e.GetMol()
        return mols

    @staticmethod
    def _strip_map(mols: list[Chem]) -> list[Chem]:
        """Removes atoms map numbers.

        Parameters
        ---------
        mols : list[Chem]
            List of molecules to clean up.

        Returns
        -------
        mols : list[Chem]
            List of molecules stripped of atom map numbers.
        """
        for m in mols:
            [a.ClearProp('molAtomMapNumber')
             for a in m.GetAtoms() if a.HasProp('molAtomMapNumber')]
        return mols


if __name__ == "__main__":
    import doctest
    doctest.testmod()
