# Mapper

**Mapper** is an implementation of approximate structure-matching algorithm
proposed by Lynch and Willett [1]. It applies the Morgan algorithm [2] with
different stopping rule *simultaneously* to both reacting molecules to detect
inter-molecular equivalences. 

The authors of Ref. [1] are a bit vague as far as assigning initial extended
connectivity  indices is concerned. Thus, besides original Morgan's approach,
two alternative methods [3,4] were implemented.

## Requirements

**Mapper** depends on:

* [Python](https://www.python.org/) (version 3.11+)
* [numpy](https://numpy.org/)
* [RDKit](http://www.rdkit.org) toolkit

## Installation

Create a virtual environment with

```
python3 -m venv <ENV_DIR>
```

where ``<ENV_DIR>`` is any directory, for example ``$HOME/venvs/mapper``.

Clone the package repository with

```
git clone https://github.com/mxk62/Mapper.git
```

Switch to the package directory and activate the virtual environment you created
earlier

```
cd Mapper
source $HOME/venvs/mapper
```

and install the package using `pip`

```
pip install --user .
```

If you want to be able to generate package documentation, you need to install
additional (optional) dependencies with


```
pip install --user .[docs]
```

instead. Then you can generate the documentation simply with

```
cd docs
make html
```

If building documentation finished successfully, open `_build/html/index.html`
in your web browser.

### Note

A fair warning though. At the moment the documentation is just an API
reference. Having said that this should provide some insight how to use the
package.

## Misc

Other files you can find in the package directory are:

* 	`demo.py`: a (very) simple CLI for demonstration purposes.
* 	`rxns.smi`: example reactions from Ref. [1] encoded in SMILES format.


## References

1.	MF Lynch and P Willett, *The Automatic Detection of Chemical Reaction
	Sites*, J Chem Inf Comput Sci **18**, 154-159 (1978)
2.	HL Morgan, , *The generation of a Unique Machine Description for Chemical
	Structures---A Technique Developed at Chemical Abstract Service*, J Chem
	Doc **5**, 107-113 (1965)
3.	CA Shelley and ME Munk, *Computer Perception of Topological Symmetry*,
	J Chem Inf Mod **17**, 110-113 (1977)
4.	K Funatsu, T Endo, N Kotera, and SI Sasaki, *Automatic Recognition of
	Reaction Site in Organic Chemnical Reactions*, Tetrahedron Computer
	Methodology **1**, 53-69 (1988)
