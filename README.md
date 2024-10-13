# RNA 3D Motif Correspondence Server (R3DMCS)

## Background
This repository contains the code for the **[R3DMCS web server](https://rna.bgsu.edu/correspondence)**<br><br>
The **RNA 3D Motif Correspondence Server** allows querying **[RNA 3D motifs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3854523/)** or specific nucleotide sets in an RNA-containing PDB structure. It then identifies corresponding nucleotides in other structures of the same molecule type from the same species or different species. The results are presented in a table based on geometric similarity, with more similar instances grouped together. An interactive heatmap with a 3D viewer visually displays the geometric differences between the motif instances. The formation of separate clusters in the heatmap may indicate variable motif geometries linked to biological function or ligand binding. R3DMCS is pronounced "Red Max".

## Documentation
Please refer to the **[Help](https://bgsu-rna.github.io/RNA-3D-correspondence/help)** page for detailed guidelines on how to create and run queries on the R3DMCS server.

Below is a brief description of the key files in the repository:
- **app.py**: The main program that manages routes in the Flask server.
- **process_input.py**: Determines the type of input provided to the server.
- **query_service.py**: Returns the complete list of nucleotides (nts) for the input query.
- **equivalence_class_service.py**: Retrieves the members of an Equivalence Class (EC) that is part of the [Representative Sets of RNA 3D structures](https://rna.bgsu.edu/rna3dhub/nrlist) for the given resolution and experimental method.
- **pairwise_service.py**: Provides pairwise annotations for equivalent nucleotide (nt) or motif instances.
- **rotation.py**: Returns the base rotation data for equivalent instances.
- **center_py**: Provides the base center data for equivalent instances.
- **discrepancy_py**: Calculates the discrepancy between equivalent nucleotide (nt) or motif instances.
- **ordering_similarity.py**: Orders equivalent nucleotide (nt) or motif instances based on similarity.
- **get_neighboring_chains.py**: Finds the neighboring chains that are within 10 Angstroms of the input nts.
- **utility.py**: Handles various utility tasks.
- **map_across_species.py**: Returns equivalent nts across Rfam family.

## Authors
- **[Sri Devan Appasamy](https://www.ebi.ac.uk/people/person/sri-devan-appasamy/)**
- **[Craig L. Zirbel](https://www.bgsu.edu/arts-and-sciences/mathematics-and-statistics/faculty-and-staff/craig-zirbel.html)**

## License
Licensed under the Apache License, Version 2.0. Please see [LICENSE](https://github.com/BGSU-RNA/RNA-3D-correspondence/blob/master/LICENSE).

## Acknowledgements
We dedicate this server to the memory of **[Neocles Leontis](https://www.bgsu.edu/arts-and-sciences/chemistry/faculty/neocles-b-leontis.html)**, our research mentor and collaborator, who inspired all of our work with RNA 3D structures.
