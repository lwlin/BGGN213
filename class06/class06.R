library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor", col="green")
points(s2.b, type="l", col="blue")
points(s3.b, type="l", col="red")

hc<- hclust(dist(rbind(s1.b,s2.b,s3.b)))
plot(hc)

#function read any type of pdb protein structure and plot on bio3d
read<- function(x){
  #4 letter PDB code
  name<- read.pdb(x)
  #trim pdb on chain A
  name.chainA<- trim.pdb(name, chain="A", elety="CA")
  #select b column of trimmed pdb file
  name.b<- name.chainA$atom$b
  #plotb3 of original structure, with secondary trimmed structrue, in line format.
  plotb3(name.b,sse=name.chainA, typ="l", ylab="Bfactor")
}
#in place of code, put your 4 letter pdb code
read("code")




#multiple proteins
#input 4ake,1ake,1e4y


#output graph with different colors








