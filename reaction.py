"""Implementation of Lynch-Willet algorithm."""

from rdkit import Chem
from chemical import Chemical


class Reaction:
    """A class representing a chemical reaction."""

    def __init__(self, smiles):
        self.reactant_smiles, self.product_smiles = smiles.split('>>')
        try:
            self.reactant = Chemical(self.reactant_smiles)
        except ValueError:
            raise ValueError("Invalid reactant SMILES")

        try:
            self.product = Chemical(self.product_smiles)
        except ValueError:
            raise ValueError("Invalid product SMILES")

    def find_core(self):
        """Extracts a reaction core.

        Function uses the Lynch-Willet algorithm to detect the reaction center.
        """

        while True:
            # Assign initial EC values to the reactant and the product.
            reactant_ecs = self.reactant.calc_init_ecs(index_type='funatsu')
            self.reactant.update_ecs(reactant_ecs)

            product_ecs = self.product.calc_init_ecs(index_type='funatsu')
            self.product.update_ecs(product_ecs)

            # If there is no common indices, exit as there is nothing to do.
            if not set(reactant_ecs).intersection(product_ecs):
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
            # which $EC_{r_{i}}^{n} = EC_{p_{j}}^{n}$.
            test_ec_order = self.reactant.ec_order
            while True:
                test_reactant_ecs = self.reactant.calc_next_ecs()
                test_product_ecs = self.product.calc_next_ecs()

                test_ec_order += 1
                common_ecs = set(test_reactant_ecs).intersection(test_product_ecs)

                if not common_ecs or test_ec_order == size_limit:
                    break

                self.reactant.update_ecs(test_reactant_ecs)
                self.product.update_ecs(test_product_ecs)

            # Find out EC-based maximal common substructure (EC-MCS).
            #
            # (1) Create maps between EC and atom indices for both product and
            # reactant.
            reactant_map = self.make_ec_map(self.reactant.ec_indices)
            product_map = self.make_ec_map(self.product.ec_indices)

            # (2) For all atoms sharing the same EC indices, find atoms
            # belonging to a circular substructure of the match radius
            # $r = n - 2$ (where $n$ is the order of last iteration),
            # centered on those atoms, i.e. EC-MCS.
            rad = test_ec_order - 2
            reactant_ec_mcs = set([])
            product_ec_mcs = set([])
            for ec in set(reactant_map).intersection(product_map):
                for idx in reactant_map[ec]:
                    reactant_ec_mcs.update(self.reactant.find_ec_mcs(idx, rad))
                for idx in product_map[ec]:
                    product_ec_mcs.update(self.product.find_ec_mcs(idx, rad))

            # (3) Delete all atoms in EC-MCS both from the reactant and
            # the product.
            self.reactant.remove_atoms(reactant_ec_mcs)
            self.product.remove_atoms(product_ec_mcs)

            # (4) Clear EC indices on the remaining atoms.
            self.reactant.clear_ecs()
            self.product.clear_ecs()

        return '>>'.join([Chem.MolToSmiles(self.reactant.mol),
                          Chem.MolToSmiles(self.product.mol)])

    def make_ec_map(self, ecs):
        """Returns a map between atoms ECs and their indices.

        Creates and returns a dictionary which keys are the EC indices and
        corresponding values are lists of atom indices with a given EC.
        """
        ec_map = {}
        for idx, ec in enumerate(ecs):
            ec_map.setdefault(ec, []).append(idx)
        return ec_map


if __name__ == '__main__':
    smarts = 'CC(=O)CC(C)C(CC#N)C(=O)N>>CC(=O)CC(C)C(CC#N)C#N'
    rxn = Reaction(smarts)
    # Should print something like C(=O)N>>C#N
    print rxn.find_core()
