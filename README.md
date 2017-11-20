# Swiss2GO
Basic FASTA to GO annotation app

Swiss2GO obtains [Gene Ontology](http://geneontology.org/) (GO) annotations for FASTA files of proteins. Annotations are important because they allow scientists to predict functions of novel proteins based on similarly structured, known proteins. 

Currently, Swiss2GO:
1. Ingests a FASTA file selected by the user using a file chooser GUI.
2. Queries the proteins of that file against the [SwissProt](http://www.uniprot.org/) database using the [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
  *BLAST usses modified dynamic programing to align the proteins of the FASTA file to the Swissprot database. 
3. Parses the output for Uniprot IDs, then obtains GO annotations.
4. Reports these annotations along with e-values of alignments back to the user in the app and in a CSV file. 

I built this app using [Tkinter](https://docs.python.org/2/library/tkinter.html) for graphics, [Biopython](http://biopython.org/) and [Requests](http://docs.python-requests.org/en/master/) for API calling. 

In the future, I'd like to expand this into an open source alternative to Blast2GO. *Side note* I'll probably change the name, I honestly had no idea Blast2GO was a thing when I started to build this and came up with a name! 

Check out a picture of Swiss2GO in use: 
![Swiss2GO](https://github.com/bkompa/bkompa.github.io/raw/master/images/Swiss2GO.png)
