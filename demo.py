import os
from reaction import Reaction


while True:
    os.system('clear')
    smiles = raw_input('Enter a reaction SMILES: ')
    if smiles == '':
        break
    print 'Processing...'
    rxn = Reaction(smiles)
    core = rxn.find_core()
    print
    print 'Reaction core:', core
    raw_input('Press <Enter> to continue.')