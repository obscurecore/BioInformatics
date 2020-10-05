# "Search for the longest ORF"

Project uses extensions:
 * `NumPy` to install run `pip install numpy`,
 * `biopython` to install run `pip install biopython` 
 where `pip` is package-management system 

Open  Reading  Frame  (ORF)  is  a  triplet  nucleotide  sequence  that  is  read  as  a  codon  that  determines  amino  acids,  one  DNA  strand  has  three  possible  reading  frames.  The  length  of  ORF  can  indicate the coding region of the candidate protein in the DNA sequence.

`ATG`  is  the  most  common  start  codon  in  DNA.

Termination codon `TAA, TAG, TGA on DNA`.


## Algorithm

 The  divide-and-conquer  algorithm  describes  problems  into  sub-problems,  recursively  sub-problems,  and  then  combines  each  solution  to  solve  the  original  problem.  However,  dynamic  programming  applies  when  sub-problems  overlap  (conditions  where  sub-problems  are  divided into the same sub-problems). Dynamic programming algorithms solve each sub-problem only once and then store the answer in a place, by avoiding the work of recomputing the answer each time completing each sub-problem [7].