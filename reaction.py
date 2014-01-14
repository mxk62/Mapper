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
        old_reactant_ecs = []
        old_product_ecs = []
        while True:
            # Assign initial EC values to the reactant and the product.
            reactant_ecs = self.reactant.calc_init_ecs(index_type='shelley')
            self.reactant.update_ecs(reactant_ecs)

            product_ecs = self.product.calc_init_ecs(index_type='shelley')
            self.product.update_ecs(product_ecs)

            # Check if EC indices changed from the last iteration. If not,
            # it means no atoms could be removed, terminate immediately.
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
            # which $EC_{r_{i}}^{n} = EC_{p_{j}}^{n}$.
            test_ec_order = self.reactant.ec_order
            while True:
                test_reactant_ecs = self.reactant.calc_next_ecs()
                test_product_ecs = self.product.calc_next_ecs()

                test_ec_order += 1
                common_ecs = set(test_reactant_ecs) & set(test_product_ecs)

                if not common_ecs or test_ec_order == size_limit:
                    break

                print test_reactant_ecs
                print test_product_ecs
                print common_ecs

                self.reactant.update_ecs(test_reactant_ecs)
                self.product.update_ecs(test_product_ecs)

            print 'Iteration stopped at {0}'.format(test_ec_order)
            print

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

            # According to Lynch and Willett, allowing match radius smaller
            # than 2 bonds significantly increases the number of incorrect
            # results. Thus, we are skipping finding the EC-MCS and
            # terminating the procedure in such a case.
            if rad < 2:
                break

            reactant_ec_mcs = set([])
            product_ec_mcs = set([])
            for ec in set(reactant_map) & set(product_map):

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

                # If there are multiple equivalents, make an arbitrary
                # assignment for each member and find relevant EC-MCSes.
                #
                # Here, we are assigning successive reactant atoms to
                # successive product atoms until all possible assignments
                # are made, i.e whichever list ends first.
                #threshold = min(len(reactant_map[ec]), len(product_map[ec]))
                #for idx in reactant_map[ec][:threshold]:
                #    reactant_ec_mcs.update(self.reactant.find_ec_mcs(idx,
                # rad))
                #for idx in product_map[ec][:threshold]:

                # Proceed to next EC index, if mapping is ambiguous.
                if len(reactant_map[ec]) != len(product_map[ec]):
                    continue

                for idx in reactant_map[ec]:
                    reactant_ec_mcs.update(self.reactant.find_ec_mcs(idx, rad))
                for idx in product_map[ec]:
                    product_ec_mcs.update(self.product.find_ec_mcs(idx, rad))

            print "Match radius:", rad
            print "R atoms slated for removal:", sorted(reactant_ec_mcs)
            print "P atoms slated for removal:", sorted(product_ec_mcs)

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

    @staticmethod
    def make_ec_map(ecs):
        """Returns a map between atoms ECs and their indices.

        Creates and returns a dictionary which keys are the EC indices and
        corresponding values are lists of atom indices with a given EC.
        """
        ec_map = {}
        for idx, ec in enumerate(ecs):
            ec_map.setdefault(ec, []).append(idx)
        return ec_map


if __name__ == '__main__':
    # Should print CC(=O)N>>CC#N.
    smarts = 'CC(=O)CC(C)C(CC#N)C(=O)N>>CC(=O)CC(C)C(CC#N)C#N'
    rxn = Reaction(smarts)
    print rxn.find_core()
