"""Class representing molecular pattern.
"""


import itertools
import re
from rdkit import Chem
from rdkit.Chem import MCS
from rdkit.Chem import SanitizeFlags
import mapper as mp


class Pattern:
    """Represents a molecular pattern."""

    def __init__(self, smarts):
        """Initializes a pattern.

        Parameters
        ----------
        smarts : string
            Pattern encoded using Daylight SMARTS notation.

        Raises
        ------
        ValueError
            This exception is raised if the input SMARTS is invalid.

        Examples
        --------
        >>> p = mp.Pattern('[C:1][O:2]')

        Pattern can be also view as an *objects* (molecular fragments), thus

        >>> all(m.HasSubstruct(p) for m, p in zip(p.fragments, p.templates))
        True
        """
        self.smarts = smarts
        self.smiles = self._strip_map(smarts)

        self.templates = [Chem.MolFromSmarts(s) for s in smarts.split('.')]
        if None in self.templates:
            raise ValueError('erroneous core SMARTS')

        self.fragments = [Chem.MolFromSmiles(s, sanitize=False)
                          for s in smarts.split('.')]
        opts = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        for m in self.fragments:
            Chem.SanitizeMol(m, sanitizeOps=opts)

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.smiles)

    def __str__(self):
        return self.smiles

    def does_contain(self, other):
        """Returns True if a pattern contains the other one.

        A pattern is considered to contain another one if the latter is the
        substructure of the first one.

        Parameters
        ----------
        other : Pattern
            Another pattern

        Examples
        --------
        Core ``[C:1][O:2]`` contain core ``[C:2]`` as the latter core's is
        a substructure of ``[C:1][O:2]``.

        >>> patt = mp.Pattern('[C:1][O:2]')
        >>> other = mp.Pattern('[C:2]')
        >>> patt.does_contain(other)
        True

        Conversely, it does not contain ``[S:2]`` as it is not a substructure
        of ``[C:1][O:2]``.

        >>> patt = mp.Pattern('[C:1][O:2]')
        >>> other = mp.Pattern('[S:2]')
        >>> patt.does_contain(other)
        False

        Multi-fragment pattern is considered to contain other
        multiple-fragment pattern if and only if:

        1. number of fragements agrees in both patterns;
        2. for *any* ordering of the latter pattern's fragments,
           *simultaneously* all of them are substructures of the former
           pattern's fragments.

        Thus

        >>> patt = mp.Pattern('[C:4]O.[N:3][C:1]=[O:2]')
        >>> other = mp.Core('[N:3].[C:2]O')
        >>> patt.does_contain(other)
        True

        but

        >>> patt = mp.Core('[C:4]O.[N:3][C:1]=[O:2]')
        >>> other = mp.Core('[C:1]=O')
        >>> patt.does_contain(other)
        False
        """
        if len(self.fragments) != len(other.fragments):
            return False
        for mols in itertools.permutations(self.fragments):
            if all(m.HasSubstructMatch(p)
                   for m, p in zip(mols, other.templates)):
                return True
        return False

    def find_distance(self, other):
        r"""Finds similarity distance between two molecular patterns.

        Similarity distance of two *connected* molecules :math:`m_{1}` and
        :math:`m_{2}` is a metric defined as

        .. math::
           d(m_{1}, m_{2}) = 1 -
             \frac{|\text{mcs}(m_{1}, m_{2}|}{\max(|m_{1}|, |m_{2}|)}

        where :math:`|\ldots|` denotes number of atoms in a given molecule.

        For any molecules :math:`m_{1}`, :math:`m_{2}` and :math:`m_{3}`,
        the following properties hold true:

        1. :math:`0 \leq d(m_{1}, m_{2}) \le 1`,
        2. :math:`d(m_{1}, m_{2}) = 0` if and only if :math:`m_{1}` and
           :math:`m_{2}` are identical,
        3. :math:`d(m_{1}, m_{2}) = d(m_{2}, m_{1})`,
        4. :math:`d(m_{1}, m_{3}) \leq d(m_{1}, m_{2}) + d(m_{2}, m_{3})`

        Function extends this concept to disconnected molecules :math:`M_{1} =
        \{m_{1}^{1}, \ldots, m_{N}^{1}\}` and :math:`M_{2} = \{m_{1}^{2},
        \ldots, m_{N}^{2}\}` for which

        .. math::
           D(M_{1}, M_{2}) =
             \begin{cases}
                \frac{1}{N} \sum_{i} d(m_{i}^{1}, m_{i}^{2}) &
                   \text{if } \forall i \quad d(m_{i}^{1}, m_{i}^{2}) < 1, \\
                   1 & \text{otherwise}.
             \end{cases}

        Parameters
        ----------
        other : Pattern
            Another pattern

        Returns
        -------
        distance : float
            The smallest similarity distance between fragments of molecular
            patterns.
        """
        if len(self.fragments) != len(other.fragments):
            return 1.0

        # At them moment, the total distance for disjoint patterns is an
        # arithmetic mean of distances between their fragments, unless
        # any pair has no MCS. In such a case it is set to 1.
        distances = []
        for mols in itertools.permutations(self.fragments):
            dfrag = []
            for mi, mj in zip(mols, other.fragments):
                mi_size, mj_size = len(mi.GetAtoms()), len(mj.GetAtoms())

                # Strictly speaking, checking for subgraph isomorphism,
                # and hence MCS algorithm, requires each of compared graphs
                # to have at least two nodes. Thus, cases where any of
                # the compared patterns consists of a single atom, are treated
                # separately using substructure matching.
                mcs_size = 0
                if min(mi_size, mj_size) != 1:
                    mcs = MCS.FindMCS([mi, mj], ringMatchesRingOnly=True)
                    if mcs.numAtoms != -1:
                        mcs_size = mcs.numAtoms
                else:
                    i, j = self.fragments.index(mi), other.fragments.index(mj)
                    if (mi.HasSubstructMatch(other.templates[j]) or
                            mj.HasSubstructMatch(self.templates[i])):
                        mcs_size = 1
                dfrag.append(1.0 - float(mcs_size) / max(mi_size, mj_size))
            dtot = 1.0 if any(abs(d - 1.0) < 1.0e-06 for d in dfrag) \
                else sum(dfrag) / len(mols)
            distances.append(dtot)
        return min(distances)

    @staticmethod
    def _strip_map(s):
        """Removes atoms map numbers."""
        tmp = re.sub(r'\[([^]]+?):(\d+)\]', r'[\1]', s)
        return re.sub(r'\[(B|b|C|c|N|n|O|o|P|S|s|F|Cl|Br|I)\]', r'\1', tmp)
