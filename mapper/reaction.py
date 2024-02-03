"""Implementation of Lynch-Willett algorithm."""

__all__ = ["Reaction"]

from collections.abc import Sequence

from rdkit import Chem
from mapper import Chemical


class Reaction:
    """A class representing a chemical reaction.

    Parameters
    ----------
    rxn_smiles : str
        A representation of a chemical reaction in Daylight SMILES
        notation.
    verbose : bool
        If true, extra messages will be provided. Defaults to False.
    """

    def __init__(self, rxn_smiles: str, verbose: bool = False) -> None:
        self.reactant_smiles, self.product_smiles = rxn_smiles.split('>>')
        self.reactant = Chemical(self.reactant_smiles)
        if self.reactant is None:
            raise ValueError('Invalid reactant SMILES')
        self.product = Chemical(self.product_smiles)
        if self.product is None:
            raise ValueError('Invalid product SMILES')
        self.verbose = verbose

    def find_core(self) -> str:
        """Extracts a reaction core.

        Function uses the Lynch-Willett algorithm to detect a reaction center.

        Returns
        -------
        core : str
            SMILES representing the reaction core.

        Examples
        --------
        Initial example from paper by Lynch and Willett:

        >>> smiles = 'CC(=O)CC(C)C(CC#N)C(=O)N>>CC(=O)CC(C)C(CC#N)C#N'
        >>> rxn = Reaction(smiles)
        >>> print(rxn.find_core())
        NC=O>>C#N
        """
        old_reactant_ecs: Sequence[int] = []
        old_product_ecs: Sequence[int] = []
        while True:
            # Assign initial EC values to the reactant and the product.
            reactant_ecs = self.reactant.calc_init_ecs(index_type='shelley')
            self.reactant.update_ecs(reactant_ecs)

            product_ecs = self.product.calc_init_ecs(index_type='shelley')
            self.product.update_ecs(product_ecs)

            # Check if EC indices changed from the last iteration. If not,
            # it means no atoms could be removed, thus terminate immediately.
            if (reactant_ecs == old_reactant_ecs and
                    product_ecs == old_product_ecs):
                break
            old_reactant_ecs = reactant_ecs
            old_product_ecs = product_ecs

            # If there is no common indices, exit as there is nothing to do.
            if not set(reactant_ecs) & set(product_ecs):
                break

            # Set size threshold to number of atoms in the smaller molecule.
            #
            # As an atom's EC index of $n$-th order can be treated as a "hash"
            # of a circular substructure of radius $(n - 1)$ bonds centered on
            # that atom, there is little point with going beyond the size of
            # the smaller molecule. Though further iterations will keep
            # increasing the EC indices on reactant and product, they mutual
            # relations should stay unchanged.
            size_limit = min(len(reactant_ecs), len(product_ecs))

            # Calculate higher order ECs until there are NO pairs of atoms for
            # which $EC_{r_{i}}^{n} = EC_{p_{j}}^{n}$ or EC order exceeds
            # the limit.
            test_ec_order = self.reactant.ec_order
            while True:
                test_reactant_ecs = self.reactant.calc_next_ecs()
                test_product_ecs = self.product.calc_next_ecs()
                test_ec_order += 1

                common_ecs = set(test_reactant_ecs) & set(test_product_ecs)
                if not common_ecs or test_ec_order >= size_limit:
                    break

                self.reactant.update_ecs(test_reactant_ecs)
                self.product.update_ecs(test_product_ecs)

            # Find out EC-based maximal common substructure (EC-MCS).
            #
            # (1) Create maps between EC and atom indices for both product and
            # reactant and find EC indices common to both the reactant and
            # the product.
            reactant_map = self.make_ec_map(self.reactant.ec_indices)
            product_map = self.make_ec_map(self.product.ec_indices)
            common_ecs = set(reactant_map) & set(product_map)

            # (2) For all atoms sharing the same EC indices, find atoms
            # belonging to a circular substructure of the match radius
            # $r = n - 2$ (where $n$ is the order of last iteration),
            # centered on those atoms, i.e. EC-MCS.
            rad = test_ec_order - 2

            if self.verbose:
                print 'Mappings (match radius: {0}):'.format(test_ec_order - 2)
                self.show_mappings(common_ecs)

            # According to Lynch and Willett, allowing match radius smaller
            # than 2 bonds significantly increases the number of incorrect
            # results. Thus, we are skipping finding the EC-MCS and terminating
            # the procedure in such a case.
            if rad < 2:
                break

            reactant_ec_mcs = set([])
            product_ec_mcs = set([])
            for ec in common_ecs:
                # If a given EC value corresponds to multiple atoms,
                # check their equivalence, i.e. all reactant atoms are
                # equivalent one to each other as well as to the product set.
                # Move to next EC value if atoms are NOT equivalent.
                reactant_classes = \
                    self.reactant.get_atom_classes(reactant_map[ec])
                product_classes = \
                    self.product.get_atom_classes(product_map[ec])
                if (len(reactant_classes) != 1 or len(product_classes) != 1 or
                        reactant_classes != product_classes):
                    continue

                # Proceed to next EC index, if mapping is ambiguous.
                if len(reactant_map[ec]) != len(product_map[ec]):
                    continue

                for idx in reactant_map[ec]:
                    reactant_ec_mcs.update(self.reactant.find_ec_mcs(idx, rad))
                for idx in product_map[ec]:
                    product_ec_mcs.update(self.product.find_ec_mcs(idx, rad))

            # (3) Delete all atoms in EC-MCS both from the reactant and
            # the product.
            self.reactant.remove_atoms(reactant_ec_mcs)
            self.product.remove_atoms(product_ec_mcs)

            # (4) Clear EC indices on the remaining atoms before next
            # iteration as they are no longer valid (atoms bond patterns may
            # have changed due to removal).
            self.reactant.clear_ecs()
            self.product.clear_ecs()

        return '>>'.join([Chem.MolToSmiles(self.reactant.mol),
                          Chem.MolToSmiles(self.product.mol)])

    def show_mappings(self, common_ecs: set[int]) -> None:
        """Prints out all atom mappings.

        Parameters
        ----------
        common_ecs : Sequence[int]
            Common EC indices.
        """
        for ec in common_ecs:
            lhs = ', '.join(['{0}:{1}'.format(*pair)
                             for pair in self._map_symbols(self.reactant, ec)])
            rhs = ', '.join(['{0}:{1}'.format(*pair)
                             for pair in self._map_symbols(self.product, ec)])
            print(ec, ': ', lhs, '<-->', rhs)

    @staticmethod
    def _map_symbols(chem: Chemical, ec: int) -> list[tuple[str, int]]:
        """Returns pairs of atom symbols and their molecular indices.

        Parameters
        ----------
        chem : Chem
            A molecule to create mapping for.

        Returns
        -------
        pairs : list[tuple[str, int]]
            List of atoms symbols paired with their molecular indices.
        """
        indices = [chem.mol.GetAtomWithIdx(i).GetProp('initIdx')
                   for i, v in enumerate(chem.ec_indices) if v == ec]
        symbols = [chem.mol.GetAtomWithIdx(i).GetSymbol()
                   for i, v in enumerate(chem.ec_indices) if v == ec]
        return list(zip(symbols, indices))

    @staticmethod
    def make_ec_map(ecs: Sequence[int]) -> dict[int, list[int]]:
        """Create a map between atoms' EC indices and their indices.

        Creates and returns a dictionary which keys are the EC indices and
        corresponding values are lists of atom indices with a given EC.

        Parameters
        ----------
        ecs : Sequence[int]
            List of EC indices.

        Returns
        -------
        ec_map : dict[int, list[int]]
            The mapping between ECs and atom indices.
        """
        ec_map = {}
        for idx, ec in enumerate(ecs):
            ec_map.setdefault(ec, []).append(idx)
        return ec_map


if __name__ == '__main__':
    # Initial example from paper by Lynch and Willett. Output should read
    #     C(=O)N>>C#N
    smiles = 'CC(=O)CC(C)C(CC#N)C(=O)N>>CC(=O)CC(C)C(CC#N)C#N'
    rxn = Reaction(smiles)
    print rxn.find_core()
