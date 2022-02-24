Atomic Data readme
These files come from Gray (3nd Ed) Appendix D1 and D2.  Ryan Terrien is responsible for scanning in the partition functions. 

There are three partition function files, to facilitate however you want to read them in.  One contains the element names, and another with the partition funcitons as tabulated in Gray.  Compare these files to Appendix D2 to confirm that you understand them.  A third file combines the information, but is somewhat trickier to parse because it contains both strings and floats.

Columns for the ionization table are:

1) Atomic number

2) Element

3) Average atomic weight

4-6) 1st, 2nd, and 3rd ionization energies in eV.

This file is most easily interpreted as being fixed width, meaning that you cannot trust whitespace to delimit the columns (some entries are blank) but you can trust the columns to be aligned.  

I also include nist_ioniz.txt, which contains much more precise ionization energies, but only for first ionization.  This is useful for getting hydrogen right.