---
layout: default
title: Help Page
---

# Webserver Help Documentation

## Table of Contents
- [Input Page](#input_page)
- [URL Access](#url_access)
---

## Input page {#input_page}
This section provides an overview of the webserver and its features.

### Step 1 - Select nucleotides
First specify the nucleotides in one 3D structure that you would like to work with.  There are different ways to describe your selection: residue number, loop id, and unit id. A short description for each selection type is given below. 

#### Residue number
For the most part, each residue in a 3D structure file can be identified by the PDB identifier (like 5J7L), the author-assigned chain identifier (like AA), and residue numbers, also called nucleotide numbers (like 1405 and 1496).  Note that the chain identifier is case sensitive, but the PDB identifier is not.

##### Individual residues
To retrieve the specific nucleotides mentioned above, one would type 1405,1496 in the Selection box, 5J7L in the PDB ID box, and AA in the Chain ID box.  The general format for entering individual nucleotide numbers is “**number1,number2,number3**" and repeat for as many individual positions as are needed.  Individual residue numbers are separated by commas.  These nucleotides are one of the closing basepairs in the decoding loop in the SSU of E. coli.  [Link to example of individual residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405,1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).  Note that in that example, some 3D structures model C1496 in syn, creating a separate cluster in the heat map.

##### Single range of residues
A range of residue numbers can be provided, separating the lower and upper number with a colon character.  The format for entering a single range of nucleotide numbers is “**start_position:end_position**”.  For example, to specify the lower-numbered strand of the decoding loop in the E. coli SSU from 5J7L chain AA, one would type 1405:1409 in the Selection box.  [Link to example of range of residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).

##### Multiple ranges of residues
This option consists of entering multiple single ranges of nucleotide numbers separated by commas. This is especially helpful when the motif is an internal loop, a junction loop, or a long-range interaction motif.  The general format is “start1:end1,start2:end2” and repeat for as many ranges as are needed.  For example, to get both strands of the E. coli SSU decoding loop, one would type 1405:1409,1491:1496.  [Link to example of multiple ranges of residue numbers](http://rna.bgsu.edu/correspondence/comparison?selection=1405:1409,1491:1496&pdb=5J7L&chain=AA&exp_method=all&resolution=3.0&scope=EC&input_form=True).

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
