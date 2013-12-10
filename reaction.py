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

        # Assign initial EC values to the reactant and the product.
        next_reactant_ecs = self.reactant.init_ecs()
        next_product_ecs = self.product.init_ecs()

        # Calculate higher order ECs until there are NO pairs of atoms for
        # which $EC_{r_{i}}^{n} = EC_{p_{j}}^{n}$.
        while set(next_reactant_ecs).intersection(next_product_ecs):
            # Save current sets of ECs.
            prev_reactant_ecs = list(next_reactant_ecs)
            prev_product_ecs = list(next_product_ecs)

            # Calculate test set of ECs for reactant and product.
            next_reactant_ecs = self.reactant.increment_ecs()
            next_product_ecs = self.product.increment_ecs()

        # Convert reactant and product to RDKit's editable molecules.
        editable_reactant = Chem.EditableMol(self.reactant.mol)
        editable_product = Chem.EditableMol(self.product.mol)

        # Delete all atoms in EC-MCS from the reactant and the product.
        reactant_map = self.make_ec_map(prev_reactant_ecs)
        product_map = self.make_ec_map(prev_product_ecs)
        for ec in reactant_map.intersection(product_map):
            [editable_reactant.RemoveAtom(i) for i in reactant_map[ec]]
            [editable_product.RemoveAtom(i) for i in product_map[ec]]



        self.reactant.remove_atoms()
        self.product.remove_atoms()

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
    smarts = 'CC(=C)CC(C)C(CC#N)C(=C)N>>CC(=C)CC(C)C(CC#N)C#N'
    rxn = Reaction(smarts)
    print rxn.find_core()
