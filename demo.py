import os

from mapper import Reaction


while True:
    os.system('clear')
    smiles = input('Enter a reaction SMILES: ')
    if smiles == '':
        break
    print('Processing...')
    rxn = Reaction(smiles, verbose=True)
    core = rxn.find_core()
    print()
    print('Reaction core:', core)
    input('Press <Enter> to continue.')
