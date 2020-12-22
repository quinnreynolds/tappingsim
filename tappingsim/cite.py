def journal_plaintext():
    return ('Olsen, J.E., Reynolds, Q.G. (2020). Mathematical Modeling of '
            'Furnace Drainage While Tapping Slag and Metal Through a Single '
            'Tap-Hole. Metall Mater Trans B 51, 1750–1759 (2020). '
            'https://doi.org/10.1007/s11663-020-01873-1')

def journal_bibtex():
    return ('@Article{Olsen2020, author={Olsen, Jan Erik and Reynolds, Quinn '
            'Gareth}, title={Mathematical Modeling of Furnace Drainage While '
            'Tapping Slag and Metal Through a Single Tap-Hole}, '
            'journal={Metallurgical and Materials Transactions B}, '
            'year={2020}, month={Aug}, day={01}, volume={51}, number={4}, '
            'pages={1750-1759}, issn={1543-1916}, '
            'doi={10.1007/s11663-020-01873-1}, '
            'url={https://doi.org/10.1007/s11663-020-01873-1}}')

def software_plaintext():
    raise NotImplementedError()

def software_bibtex():
    raise NotImplementedError()
