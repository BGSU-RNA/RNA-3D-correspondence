---
layout: default
title: Help Page
---

# Webserver Help Documentation

## Table of Contents
- [Input Page](#input_page)
- [URL Access](#url_access)
- [Output page and examples](#output_page)


---

## Input page {#input_page}
This section provides an overview on how to fill up the input page.

### Step 1 - Select nucleotides
First specify the nucleotides in one 3D structure that you would like to work with.  There are different ways to describe your selection: residue number, loop id, and unit id. A short description for each selection type is given below. 

#### Residue number
For the most part, each residue in a 3D structure file can be identified by the PDB identifier (like 5J7L), the author-assigned chain identifier (like AA), and residue numbers, also called nucleotide numbers (like 1405 and 1496).  Note that the chain identifier is case sensitive, but the PDB identifier is not.

##### Individual residues
To retrieve the specific nucleotides mentioned above, one would type 1405,1496 in the Selection box, 5J7L in the PDB ID box, and AA in the Chain ID box.  The general format for entering individual nucleotide numbers is “**number1,number2,number3**" and repeat for as many individual positions as are needed.  Individual residue numbers are separated by commas.  These nucleotides are one of the closing basepairs in the decoding loop in the SSU of E. coli.  [Link to example of individual residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).  Note that in that example, some 3D structures model C1496 in syn, creating a separate cluster in the heat map.

##### Single range of residues
A range of residue numbers can be provided, separating the lower and upper number with a colon character.  The format for entering a single range of nucleotide numbers is “**start_position:end_position**”.  For example, to specify the lower-numbered strand of the decoding loop in the E. coli SSU from 5J7L chain AA, one would type 1405:1409 in the Selection box.  [Link to example of range of residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).

##### Multiple ranges of residues
This option consists of entering multiple single ranges of nucleotide numbers separated by commas. This is especially helpful when the motif is an internal loop, a junction loop, or a long-range interaction motif.  The general format is “**start1:end1,start2:end2**” and repeat for as many ranges as are needed.  For example, to get both strands of the E. coli SSU decoding loop, one would type 1405:1409,1491:1496.  [Link to example of multiple ranges of residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).

#### Loop ID
Each week, the BGSU data processing pipeline extracts hairpin, internal, and junction loops from RNA-containing 3D structure files using the FR3D software. Once the loops are extracted, we label them with unique and stable identifiers. These “loop ids" contain the following three fields, separated by underscores:
- Field 1: Loop type prefix: “HL” for hairpin loops, “IL” for internal loops, “J3” for three-way junctions
- Field 2: PDB ID
- Field 3: A sequentially assigned, three digit character

Users can view the loop_ids for a particular RNA structure by exploring the RNA Structure Atlas pages.  See, for example, the [page for 5TBW](http://rna.bgsu.edu/rna3dhub/pdb/5TBW), hairpin loop [HL_5TBW_007](http://rna.bgsu.edu/rna3dhub/loops/view/HL_5TBW_007), internal loop [IL_5TBW_019](http://rna.bgsu.edu/rna3dhub/loops/view/IL_5TBW_019), and 3-way junction loop [J3_5TBW_003](http://rna.bgsu.edu/rna3dhub/loops/view/J3_5TBW_003).  Loop ids are also used in the [RNA 3D Motif Atlas](http://rna.bgsu.edu/rna3dhub/motifs).  For example, motif group [IL_29549.7](http://rna.bgsu.edu/rna3dhub/motif/view/IL_29549.7) contains 35 instances of the kink-turn internal loop motif.  Each of the 35 loop ids there could be used as the starting point to explore variation in the kink turn geometry across different experimental structures. 

To specify a R3DMCS query using a loop id, type the loop id in the Selection box, and leave the PDB id and Chain id boxes empty.  [Link to example of using loop id to specify a query](http://rna.bgsu.edu/correspondence/comparison?selection=J3_5TBW_003&exp_method=all&resolution=3.0&scope=EC&input_form=True).

#### Unit ID
Unit ids allow full ability to specify a particular residue.  Unit ids are described on the [unit_id page](https://www.bgsu.edu/research/rna/help/rna-3d-hub-help/unit-ids.html).  For example, 5J7L|1|AA|G|1405 and 5J7L|1|AA|C|1496 are unit ids for one of the flanking pairs in the E. coli SSU decoding loop.  Each unit id tells the PDB id, the model number, the chain, the sequence of the residue, and the residue number, and can optionally also tell the insertion code, alternate id, and symmetry operator. [Link to example of using unit ids to specify a query](http://rna.bgsu.edu/correspondence/comparison?selection=8GLP%7C1%7CL5%7CG%7C1561,8GLP%7C1%7CL5%7CG%7C1562,8GLP%7C1%7CL5%7CA%7C1563,8GLP%7C1%7CL5%7CA%7C1564,8GLP%7C1%7CL5%7CA%7C1565,8GLP%7C1%7CL5%7CC%7C1566&input_form=True).

#### Special cases
In NMR structures, model number is sometimes used.  The residue number option above assumes that model 1 is desired, and provides no way to specify a different model.  Loop ids can come from models other than 1.  Unit ids can be used to specify PDB id, model number, chain id, residue sequence, residue number, alternate id, and symmetry operator. [Link to example of using unit ids to specify an NMR model](http://rna.bgsu.edu/correspondence/comparison?selection=2I7Z%7C4%7CA%7CU%7C5,2I7Z%7C4%7CA%7CA%7C6,2I7Z%7C4%7CA%7CU%7C7,2I7Z%7C4%7CA%7CG%7C8,2I7Z%7C4%7CA%7CC%7C35,2I7Z%7C4%7CA%7CU%7C36,2I7Z%7C4%7CA%7CA%7C37,2I7Z%7C4%7CA%7CA%7C38,2I7Z%7C4%7CA%7CG%7C39&exp_method=nmr&scope=EC&input_form=True).

### Step 2 - Correspondence Scope
Once you specify the nucleotides of interest from a particular structure, R3DMCS can retrieve nucleotides for comparison in one of two ways:

#### Corresponding nucleotides from the same molecule in the same organism
R3DMCS uses the Equivalence Classes in the [BGSU RNA Representative Sets](http://rna.bgsu.edu/rna3dhub/nrlist) to identify the 3D structures of the same molecule from the same organism.  For example, the [E. coli SSU equivalence class NR_3.0_56726.118](http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_3.0_56726.118) collects all 151 small subunit ribosomal RNA chains from different 3D structures solved at 3.0Å or better resolution, as of April 3, 2024.  Sequence-based alignments between these chains are computed each week for newly added chains, and R3DMCS uses those alignments to find corresponding nucleotides.

#### Homologous nucleotides from the same molecule in different organisms
R3DMCS uses multiple sequence alignments produced by the BGSU RNA pipeline using [Infernal](http://eddylab.org/infernal/) and based on [Rfam covariance models](https://docs.rfam.org/en/latest/glossary.html#:~:text=Clustal%20web%20server.-,Covariance%20model%20(CM),in%20RNA%20or%20DNA%20sequences.), for all PDB chain sequences in each Rfam family.  The Rfam covariance models are downloaded with each Rfam release.  Infernal and the covariance models are used each week to map PDB chains to Rfam families.  For example, the E. coli SSU chain AA from PDB structure 5J7L maps to Rfam family RF00177.  To see all of the other chains that map to Rfam family RF00177, one can go to a [Representative Set release page](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.332/3.0A) and type RF00177 in the Filter box.  As of the April 24, 2024 release at 3.0Å resolution threshold, this results in 21 equivalence classes.  Of these, 15 correspond to bacterial species, 5 to mitochondrial ribosomes, and 1 to a chloroplast ribosome.  The number of chains in each equivalence class varies from 1 to 206.  To avoid over-representation from any one equivalence class, the following input parameter can be used when retrieving over different organisms:

##### Depth
This parameter tells the maximum number of chains in each equivalence class to use for the retrieval across different organisms.  The chains are prioritized in the order that they are listed on the Representative Set page, which is according to six structure quality factors, which are briefly described on the [Representative Set home page](http://rna.bgsu.edu/rna3dhub/nrlist).

### Step 3 - Resolution threshold
The BGSU RNA site organizes equivalence classes (chains of the same molecule from the same species) by resolution, keeping all structures solved up to a threshold value.  Setting the resolution threshold makes it possible to look only at the highest resolution structures.  The query will retrieve corresponding motif instances from all chains in the equivalence class with resolution up to the resolution threshold.  The default resolution threshold is 3 Ångstroms.  Note that NMR structures have no reported resolution, and can only be retrieved by setting the resolution threshold to "all".

### Step 4 - Experimental technique
Choose the desired experimental technique. Users can choose X-ray diffraction, cryo-EM microscopy, or “all” which also includes NMR and other techniques besides X-ray and cryo-EM.

### Step 5 - PDB identifiers to exclude
Occasionally a PDB entry will have a very different 3D structure in a specific region, due to a variety of factors such as the presence of a modified nucleotide, a bound ligand, or a modeling feature that sets it far apart from other structures.  If it is desired to exclude those instances in order to focus on ones of more direct interest, it is possible to give a list of PDB identifiers to exclude.  For example, in the [example of individual residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True) from above, the structures 4V9O and 4V9P model C1496 in syn, which sets them apart from the other instances.  Also, structure 8GHU has the non-standard residue ZIV at position 1405, which has two bases in one nucleotide and so is not compared geometrically to the other instances, which is shown with gray cells in the heat map.  [Listing these three structures under PDB identifiers to exclude](http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&exclude=4V9O,4V9P,8GHU&input_form=True) removes these instances, and the heat map coloring adjusts accordingly, making it easier to see the high degree of similarity between the remaining instances.

### Step 6 - Submit
After pressing the Submit button, the R3DMCS server will retrieve data for all corresponding nucleotides, calculate all-against-all discrepancies, order the instances by similarity, and render the page in HTML.  This process typically takes 5 to 15 seconds, depending on the number of instances returned.  However, a larger number of instances will take longer; for reference, R3DMCS took almost 6 minutes for a tRNA hairpin retrieval that resulted in 718 instances.  It seems that that was an unnecessarily broad query.

## URL access {#url_access}
R3DMCS works by responding to a carefully created URL.  The only function of the input page described above is to create the URL, and then R3DMCS processes the URL.  The URL can be created separately from the input page using the settings described in this section.

### Base URL
The base URL is:
- [http://rna.bgsu.edu/correspondence/comparison](http://rna.bgsu.edu/correspondence/comparison)

### Selection of nucleotides
It is required to specify a selection, using one of the methods listed above.  For example, to specific a hairpin loop in PDB structure 5TBW:
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007)
Default settings in the query above are to retrieve within the equivalence class, resolution threshold 3.0Å, all experimental techniques, and no excluded PDB ids.

To use residue numbers, also specify the PDB id and the chain id, and then specify ranges as explained in the examples above.  For individual nucleotides:
- [http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&resolution=3.0](http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&resolution=3.0)

For a single range of nucleotides:
- [http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409&pdb=5J7L&chain=AA&resolution=3.0](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409&pdb=5J7L&chain=AA&resolution=3.0)

For multiple ranges of nucleotides:
- [http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496&pdb=5J7L&chain=AA&resolution=3.0
](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496&pdb=5J7L&chain=AA&resolution=3.0)

For multiple ranges of nucleotides plus an individual nucleotide:
- [http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496,530&pdb=5J7L&chain=AA&resolution=3.0](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496,530&pdb=5J7L&chain=AA&resolution=3.0)

Note that we are including the resolution threshold in these examples, to provide smaller and faster queries.

### Scope and Depth
To specify the scope use the "scope" key and values EC for retrievals across the equivalence class of same molecule, same species, and use Rfam for retrievals across species.
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC)
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=Rfam](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=Rfam)

When using scope=Rfam, the default depth is 1, but other values can be chosen, so that R3DMCS retrieves additional instances from each species:
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=Rfam&depth=4](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=Rfam&depth=4)

### Resolution threshold
To specify the maximum resolution, use the "resolution" key and value chosen from 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 20.0, all:
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0)
Make sure that the query structure is within the resolution threshold.

### Experimental technique
To specify the experimental technique, use the exp_method key and values all, x-ray, xray, cryo-em, em, nmr:
- [http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0&exp_method=x-ray](http://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0&exp_method=x-ray)
- [https://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0&exp_method=cryo-em](https://rna.bgsu.edu/correspondence/comparison?selection=HL_5TBW_007&scope=EC&resolution=3.0&exp_method=cryo-em)
- [http://rna.bgsu.edu/correspondence/comparison?selection=IL_2I7Z_001&exp_method=nmr&scope=EC](http://rna.bgsu.edu/correspondence/comparison?selection=IL_2I7Z_001&exp_method=nmr&scope=EC)

Note that when exp_method is nmr, then the resolution will be set to all, so that NMR structures will be found and returned.

###   PDB IDs to exclude
To exclude specific PDB ids, for example because they have an unusual feature that distracts from the main data analysis, use the exclude key and list PDB ids separated by commas:
- [https://rna.bgsu.edu/correspondence/comparison?pdb=5J7L&chain=AA&selection=1405%3A1409&exp_method=all&resolution=3.0&scope=EC&exclude=8GHU,8G2U](https://rna.bgsu.edu/correspondence/comparison?pdb=5J7L&chain=AA&selection=1405%3A1409&exp_method=all&resolution=3.0&scope=EC&exclude=8GHU,8G2U)

In the URL above, we modify the earlier example of nucleotides 1405 to 1409 from chain AA of 5J7L to exclude 8GHU and 8G2U. The instance from 8GHU has the modified nucleotide ZIV in position 1405 and so no discrepancy was calculated. The instance from 8G2U had large discrepancy with nearly every other nucleotide, causing the color range to be mostly blue between all other instances. Excluding those instances makes it possible to better discern the structure of the remaining instances; the maximum discrepancy dropped from 1.26 to 0.39.

Another good example is obtained from the E. coli SSU basepair 1405 with 1496; in the query above, two instances have the C modeled in syn, which are real outliers compared to the other instances.  Excluding 4V9O and 4V9P together with 8G2U and 8GHU allows us to focus on the other instances.
- [https://rna.bgsu.edu/correspondence/comparison?pdb=5J7L&chain=AA&selection=1405%2C1496&exp_method=all&resolution=3.0&scope=EC&exclude=4V9P,4V9O,8G2U,8GHU](https://rna.bgsu.edu/correspondence/comparison?pdb=5J7L&chain=AA&selection=1405%2C1496&exp_method=all&resolution=3.0&scope=EC&exclude=4V9P,4V9O,8G2U,8GHU)

###   Pre-filled input page
To direct the user to the input page with fields already filled in, so the user can modify the fields as desired before requesting the results, use the input_form key with value true:
- [http://rna.bgsu.edu/correspondence/comparison?selection=IL_5J7L_014&exp_method=all&resolution=3.0&scope=EC&input_form=True](http://rna.bgsu.edu/correspondence/comparison?selection=IL_5J7L_014&exp_method=all&resolution=3.0&scope=EC&input_form=True)

Note that a significant fraction of the instances on this page have A279 modeled in syn rather than anti in the example above.

Links to the input page are provided on all hairpin, internal, and 3-way junction loop pages on the BGSU RNA site; see for example the [page for IL_5J7L_014](http://rna.bgsu.edu/rna3dhub/loops/view/IL_5J7L_014), which is a kink turn.

After loading the input page, it takes 2 seconds for the Submit button to change to blue and become active.  That is to slow down bots that might click links to the R3DMCS input page. 

## Output page and examples {#output_page}

