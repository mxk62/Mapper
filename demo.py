import os
from mapper.reaction import Reaction


while True:
    os.system('clear')
    smiles = raw_input('Enter a reaction SMILES: ')
    if smiles == '':
        break
    print 'Processing...'
    rxn = Reaction(smiles, verbose=True)
    core = rxn.find_core()
    print
    print 'Reaction core:', core
    raw_input('Press <Enter> to continue.')