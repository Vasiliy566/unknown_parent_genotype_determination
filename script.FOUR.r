setwd("/Users/vasiliyisaev/Desktop/Bioinformatics/Practice")
DATA=read.table("SNP.txt",stringsAsFactors=FALSE)
head(DATA)
DATA=DATA[,c(1,3,20, 28, 29, 33:38)]
head(DATA)

names(DATA)=c('uploaded_variant','allele','existing_var','given_ref','used_ref','AF',"AFR",'AMR','EAS','EUR','SAS')
uDATA=unique(DATA)
head(uDATA)
uDATA=subset(uDATA, uDATA$AF!='-')
head(uDATA)
dim(uDATA)
write.csv(uDATA, 'SNP.csv')
uDATA=read.csv('SNP.csv')
head(uDATA)
hist(uDATA$AF)
P=read.table('plink.tped')[,c(2,5:8)]; head(P)
names(P)=c('uploaded_variant','C1','C2','M1','M2')
head(P)
D=merge(P,uDATA, by='uploaded_variant')
head(D)
D=(subset(D, D$given_ref!='-'))
D=(subset(D, D$given_ref==D$used_ref))


write.csv(D,'ALL.csv')
D=read.csv('ALL.csv')
N=dim(D)[1]
complement<-function(x){
  if(x=='A'){
    return ('T');
    }
  else if(x=='T'){
    return ('A');
    }
  else if(x=='C') {
    return ('G');
    }
  else if(x=='G') {
    return ('C');
    }
  else{
    print('error in function complement with arg:' );
    print(x);
    return ('X');
  }
}

write.table(paste('marker','alleleA','alleleB','1001','1001','1001', sep=' '),
                  'EGOR.beagle.uploaded.txt', quote=F, row.names = F, col.names = F)
for(n in 1:N){
  #D[n,]
  #print('----')
  #( D[n,2]=complement(D[n,2]) )
  #print('----')
  print(n)
  #print(D[n,c(2,3,4,5,6,8)])
  if(as.character(D[n,2])!=as.character(D[n,6]) & as.character(D[n,2])!=D[n,8])
  {
    D[n,2]=complement(D[n,2]);
    D[n,3]=complement(D[n,3]);
    D[n,4]=complement(D[n,4]);
    D[n,5]=complement(D[n,5]);
    
  }
  #D[n,c(2,3,4,5,6,8)]
  #consider different cases A=reference B=alternative
  #case 1 mother heterozygous AB child AA
  #prior  frequencies of father's genotype from population
  #assuming european
  fAA=(1-D[n,14])^2; fAA;
  fAB=2*(1-D[n,14])*D[n,14]; fAB;
  fBB= (D[n,14])^2; fBB;
  PAA=0;
  PBB=0;
  PAB=0;
  #case 1 mother AA 
  if(D[n,4]==D[n,5] &
     D[n,5]==as.character(D[n,8]))
  {
    print ('CASE 1');
    #child AA
    if( D[n,2]==D[n,3] & D[n,2]==D[n,5] )
    {PBB=0;
    PAA=fAA;
    PAB=0.5*fAB;}
  #child AB
  else
  {
    PBB=fBB;
    PAB=0.5*fAB;
    PAA=0  
  }
  }
  #case 2 mother AB 
  if(D[n,4]!=D[n,5] )
  {
    print ('CASE 2');
    #child AA
    if(as.character(D[n,8])==D[n,2] & D[n,2]==D[n,3])
    {
      PBB=0;
    PAB=0.25*fAB;
    PAA=0.5*fAA;  
    }
    else if(D[n,2]!=D[n,3]){#child AB
      PBB=0.5*fBB;
      PAB=0.5*fAB;
      PAA=0.5*fAA;  
    }
    else{#child BB
      PAA=0;
      PBB=0.5*fBB;
      PAB=0.5*fAB;  
    }
  }
  #case 3 mother BB
  if(D[n,4]==D[n,5] & D[n,5]!=as.character(D[n,8]))
  {
    print ('CASE 3');
    if(D[n,2]==D[n,3]){#child BB
      PAA=0;
      PAB=fAB*0.5;
      PBB=fBB;
    }
    else{#child AB
      PAA=fAA;
      PAB=fAB*0.5;
      PBB=0;
    }
  }
   write.table(paste(D[n,]$uploaded_variant, D[n,]$given_ref, D[n,]$allele, PAA,PAB, PBB, sep=' '),'EGOR.beagle.txt', quote=F, row.names = F, col.names = F, append = T)
}

#get chromosomal positions of rs
D1$POS=as.character(D1$marker)
for(i in 1:dim(D1)[1]){
  print(i)
  RS=D1[i,1]
  I=(regexpr(",", RS, fixed = TRUE) - 1)[[1]]
  if(I>0){RS=substr(RS,1,I)}
  print(RS)
  K=annotations(snp = RS, output = 'metadata')
  D1[i,7]=as.character(paste(K[2,2],K[3,2],sep = '_'))
#annotations(snp = 'rs7903146', output = 'snpedia')
}
write.table(D1,'EGOR.pos.txt')

E=read.csv('EGOR.csv',header=T)
M=read.csv('marker_pos.csv')
EM=merge(E, M, by='marker')#[, c(7,2:6)]
#allele1 reference (major), allele2 alternative (minor)
names(EM)=c('marker','allele1', 'allele2','Ind0','Ind0','Ind0','SNP')
write.table(EM, 'EM.txt',row.names = F, quote=F)

#read reference
REF=read.table('plink.frq.strat',header=T)#A1 minor, A2 major
length(unique(REF$SNP))

REF_E=merge(REF, E, by.x='SNP', by.y='marker')
REF_E_CONS=subset(REF_E, REF_E$A2==REF_E$alleleA); dim(REF_E_CONS)
REF_E_NCONS=subset(REF_E, REF_E$A2!=REF_E$alleleA); dim(REF_E_NCONS)

dim(REF_E_CONS)
dim(REF_E_NCONS)
head(REF_E_NCONS)

#this part is not necessary I think
#REF_E_NCONS$A1=REF_E_NCONS$alleleB;
#REF_E_NCONS$A2=REF_E_NCONS$alleleA;
#REF_E_NCONS$MAF=1-REF_E_NCONS$MAF;
#head(REF_E_NCONS)
#dim(REF_E)
#REF_E=rbind(REF_E_CONS,REF_E_NCONS)
#dim(REF_E)

REF1=merge(REF_E, EM, by.x='SNP', by.y='marker')
head(REF1)
names(REF1)[1]='marker'
#for(k in 1: dim(REF_E)[1]){
#  if(REF_E[k,]$A1!=REF_E[k,]$alleleA){
##    REF_E[k,]$A1=REF_E[k,]$alleleA;
 #   REF_E[k,]$A2=REF_E[k,]$alleleB;
 #   REF_E[k,]$MAF=1-REF_E[k,]$MAF;
 # }
#}
attach(REF1)
REF1=REF1[order(marker, CLST),]
head(REF1)
#head(subset(REF_E, REF_E$A1!=REF_E$alleleA))
#id chr pos name A0_freq A1 French Han Chukchi Karitiana Papuan Sindhi Yoruba
#1_752566 1 752566 rs3094315 G A 0.166666666666667 0.0606060606060606 0.369565217391304 0.0833333333333334 0.0714285714285714 0.305555555555556 0.671428571428571

L=length(unique(REF1$SNP)); L
M=length(unique(REF1$CLST))+6;M

NEW_REF=data.frame(matrix(vector(), L, M,
                       dimnames=list(c(), 
                                     c("id","chr","pos","name","A0_freq","A1",
                                       as.character(unique(REF1$CLST))))),
                stringsAsFactors=F)
NEW_REF$id=unique(REF1$SNP)
NEW_REF[,c(1,2,4,5,6)]=unique(REF1[,c(19,2,1,4,5)])
head(NEW_REF)
NEW_REF$pos=substr(NEW_REF$id,hashPos,100)
head(NEW_REF)

names(NEW_REF)[7]

#for(col in 7:dim(NEW_REF)[2]){
#  for(row in 1:dim(NEW_REF)[1]){
# maf=subset(REF1, REF1$CLST==names(NEW_REF)[col]&REF1$marker==NEW_REF$name[row])$MAF
#NEW_REF[row, col]=maf;
#  }}
library('Rmpi')
mpi.spawn.Rslaves(nslaves=2)
mpi.spawn.Rslaves(nslaves=mpi.universe.size()-1)
#for(col in 7:7){
  for(row in 1:dim(NEW_REF)[1]){#REF1$CLST==names(NEW_REF)[col]&
    maf=subset(REF1, REF1$marker==NEW_REF$name[row])$MAF
    NEW_REF[row, 7:(dim(NEW_REF)[2])]=maf;
    #print(NEW_REF$name[row])
  }#}
#NEW_REF1=NEW_REF[,1:62]
write.table(NEW_REF, 'EGOR_REF.2.txt',quote = F, row.names = F)

#s <- strsplit(as.character(NEW_REF$id), "_")
hashPos = regexpr("_", NEW_REF$id, fixed=TRUE) + 1
dotPos = length(as.character(NEW_REF$id[1]))
finalText = substring(NEW_REF$id[1], hashPos, dotPos)
finalText


rm_between(S, ":", ";", extract=TRUE)



#AUGUST 2 2018
#new mapping of RS to positions
DAT=read.table('H_parsed.txt')
head(DAT)
names(DAT)=c('CHR','POS','DESC')
library(qdapRegex)
library(stringr)

DAT$RS=rm_between(DAT$DESC, ":", ";", extract=TRUE)
DAT$REF=rm_between(DAT$DESC, "ancestral_allele=",";",  extract=TRUE)
DAT$ALT=rm_between(DAT$DESC, "Variant_seq=",";",  extract=TRUE)
DAT$MAF=rm_between(DAT$DESC, "|","|",  extract=TRUE)

#DAT$REF=str_sub(DAT$DESC, start=-1,end = -1)
head(DAT)

#DATA_RS=read.csv('DATA_RS.csv', row.names = 1)
#dim(DATA_RS)
#VAR=read.table('var.txt')
#names(VAR)='RS'
#head(VAR)
#VAR_DATA=dim(merge(VAR, DAT, by='RS'))

#VAR_DATA$REF=rm_between(VAR_DATA$DESC, "Reference_seq=", ";", extract=TRUE)
#head(VAR_DATA)
#dim(VAR_DATA)
#VAR_DATA$marker=paste(VAR_DATA$CHR,VAR_DATA$POS, sep='_')
#VAR_FINAL=VAR_DATA[,c(1,2,3,5,6)]
#dim(VAR_FINAL)

P=read.table('plink.tped'); head(P)
names(P)=c('N0','uploaded_variant','N1','N2','C1','C2','M1','M2')
head(P)
D=as.matrix(merge(P[,c(2,5:8)],DAT[,c(1,2,4,5,6,7)], by.x='uploaded_variant', by.y='RS'))
dim(D)
head(D)


write.csv(D,'ALL.csv')
D=read.csv('ALL.csv')
N=dim(D)[1]
D$P1=D$M1
D$P2=D$M2
complement<-function(x){
  if(x=='A') return ('T');
  if(x=='T') return ('A');
  if(x=='C') return ('G');
  if(x=='G') return ('C');
}

D[is.na(D$MAF),]$MAF=0.5
D=subset(D, D$C1!=0 & D$C2!=0&D$M1!=0 &D$M2!=0&!is.na(D$POS))
N=dim(D)[1]
write.table(paste('marker','alleleA','alleleB','Ind0','Ind0','Ind0'
                  , sep=' '),
            'father.txt', quote=F, row.names = F, col.names = F)
for(n in 1:N){
  #D[n,]
 
  #consider different cases A=reference B=alternative
  #case 1 mother heterozygous AB child AA
  #prior  frequencies of father's genotype from population
  #assuming european
  maf=as.numeric(D[n,]$MAF)
  fAA=(1-maf)^2; fAA;#A -reference B alternative
  fAB=2*(1-maf)*maf; fAB;#hetero
  fBB= (maf)^2; fBB;#ancestral
  PAA=0;
  PBB=0;
  PAB=0;
  #variants
  c1=as.character(D[n,]$C1)
  c2=as.character(D[n,]$C2)
  m1=as.character(D[n,]$M1)
  m2=as.character(D[n,]$M2)
  ref=as.character(D[n,]$REF)
  alt=as.character(D[n,]$ALT)
  W=c(ref, alt )
  if(c1 %in% W & c2 %in% W & m1 %in% W & m2 %in% W){
  #case 1 mother AA 
  if(m1==m2 &m1==ref)
  {
    print ('CASE 1');
    #child AA
    if( c1==c2 &c1==ref ){
      PBB=0;
    PAA=fAA;
    PAB=0.5*fAB;}
    #child AB
    else
    {
      PBB=fBB;
      PAB=0.5*fAB;
      PAA=0  ; 
    }
  }
  #case 2 mother AB 
  if(m1!=m2)
  {
    print ('CASE 2');
    #child AA ref ref
    if(c1==ref & c1==ref){
    
      PBB=0;
      PAB=0.25*fAB;
      PAA=0.5*fAA;  
    } 
    if(c1!=c2 ){
      #child AB
      PBB=0.5*fBB;
      PAB=0.5*fAB;
      PAA=0.5*fAA;  
    
    }
    if(c1==c2& c1==alt){
      PAA=0;
      PBB=0.5*fBB;
      PAB=0.5*fAB;  
    }
  }
  #case 3 mother BB
  if(m1==m2 & m1==alt)
  {
    print ('CASE 3');
    if(c1==c2 & c1==alt){#child BB
      PAA=0;
      PAB=fAB*0.5;
      PBB=fBB;
    } else{#child AB
      PAA=fAA;
      PAB=fAB*0.5;
      PBB=0;
    }
  }
  #A ref
  #B alr
  if(PAA>PAB & PAA> PBB){ #homoref
    D[n,]$P1=ref;
    D[n,]$P2=ref;
    }
  if(PBB>PAB & PBB> PAA){ #homo alt
    D[n,]$P1=alt;
    D[n,]$P2=alt;
    }
  if(PAB>PBB& PAB>PAA){
    D[n,]$P1=alt;
    D[n,]$P2=ref;
    }
  print(PAA)
  print(PBB)
  print(PAB)
  BS=PAA+PAB+PBB;
  PAA=PAA/BS;
  PAB=PAB/BS;
  PBB=PBB/BS;
  marker=paste(D[n,]$CHR,D[n,]$POS,sep="_");
  write.table(paste(marker, ref, alt, PAA,PAB, PBB, sep=' '),'father.txt', quote=F, row.names = F, col.names = F, append = T)
  }
}

write.csv(D, "ALLPED.txt")
B=read.table('EGOR.beagle.aug2.txt', header=T)
DB=merge( B,D, by.x='marker',by.y = 'uploaded_variant')

#get mothers genotype

write.table(paste('marker','alleleA','alleleB','ind0','ind0','ind0','ind2'),
            'child.txt', quote=F, row.names = F, col.names = F)
write.table(paste('marker','alleleA','alleleB','ind0','ind0','ind0','ind2'),
            'mother.txt', quote=F, row.names = F, col.names = F)
N=dim(D)[1]; N; D$FLAG='GOOD'; head(D)
for(n in 1:N){
  #D[n,]
  #print(n)
  #print(D[n,])
    c1=as.character(D[n,]$C1)
    c2=as.character(D[n,]$C2)
    m1=as.character(D[n,]$M1)
    m2=as.character(D[n,]$M2)
    ref=as.character(D[n,]$REF)
    alt=as.character(D[n,]$ALT)
    W=c(ref, alt )
    if(!(c1 %in% W)){c1=complement(c1);D[n,]$C1=c1; print('complement')}
    if(!(c2 %in% W)){c2=complement(c2);D[n,]$C2=c2; print('complement')}
    if(!(m1 %in% W)){m1=complement(m1);D[n,]$M1=m1; print('complement')}
    if(!(m2 %in% W)){m2=complement(m2);D[n,]$M2=m2; print('complement')}
    
    if(!(c1 %in% W)|(!(c2 %in% W))| (!(m1 %in% W))|(!(m2 %in% W))){D[n,]$FLAG='BAD'}
    else{#if(D$FLAG=='GOOD'){
  #mother and child
  MAA=0.01; MAB=0.01; MBB=0.01;
  CAA=0.01; CAB=0.01; CBB=0.01;
  
  #case 1 mother AA 
  if(m1==m2 & m1==ref)
  {
   # print ('CASE 1');
    MAA=1;
  }
  #case 2 mother AB 
  if(m1!=m2 ) {
    MAB=1;
    #print ('CASE 2');
    }
  #case 3 mother BB
  if(m1==m2 & m1==alt){
    MBB=1;
   # print ('CASE 3');
    }
  
  #child
  if(c1==c2 &c1==ref) {
    #print ('CASE 1');
    CAA=1;
  }
  #case 2 mother AB 
  if(c1!=c2) {
    CAB=1;
    #print ('CASE 2');
  }
  #case 3 mother BB
  if(c1==c2 &  c1==alt ) {
    CBB=1;
    #print ('CASE 3');
  }
  if(CAA+CAB+CBB!=1.02 ||MAA+MAB+MBB!=1.02  ){print(n); break;}
  marker=paste(D[n,]$CHR,D[n,]$POS,sep='_')
  write.table(paste(marker, ref, alt, MAA,MAB, MBB, sep=' '),'mother.txt', quote=F, row.names = F, col.names = F, append = T)
  write.table(paste(marker, ref,alt, CAA,CAB, CBB, sep=' '),'child.txt', quote=F, row.names = F, col.names = F, append = T)
    }
  }
D1=subset(D, D$FLAG=='GOOD')
dim(D)
dim(D1)

subset(D, D$FLAG=='BAD')


#get chromosomal positions of rs
D1$POS=as.character(D1$marker)
for(i in 1:dim(D1)[1]){
  print(i)
  RS=D1[i,1]
  I=(regexpr(",", RS, fixed = TRUE) - 1)[[1]]
  if(I>0){RS=substr(RS,1,I)}
  print(RS)
  K=annotations(snp = RS, output = 'metadata')
  D1[i,7]=as.character(paste(K[2,2],K[3,2],sep = '_'))
  #annotations(snp = 'rs7903146', output = 'snpedia')
}
write.table(D1,'EGOR.pos.txt')

#finx D
for(n in 1:N){
  #D[n,]
  #print(n)
  #print(D[n,])
  c1=as.character(D[n,]$C1)
  c2=as.character(D[n,]$C2)
  m1=as.character(D[n,]$M1)
  m2=as.character(D[n,]$M2)
  ref=as.character(D[n,]$REF)
  alt=as.character(D[n,]$ALT)
  p1=as.character(D[n,]$P1)
  p2=as.character(D[n,]$P2)
  
  W=c(ref, alt )
  if(!(c1 %in% W)){c1=complement(c1);D[n,]$C1=c1; print('complement')}
  if(!(c2 %in% W)){c2=complement(c2);D[n,]$C2=c2; print('complement')}
  if(!(m1 %in% W)){m1=complement(m1);D[n,]$M1=m1; print('complement')}
  if(!(m2 %in% W)){m2=complement(m2);D[n,]$M2=m2; print('complement')}
  if(!(p1 %in% W)){p1=complement(p1);D[n,]$P1=p1; print('complement')}
  if(!(p2 %in% W)){p2=complement(p2);D[n,]$P2=p2; print('complement')}
  
  if(!(c1 %in% W)|(!(c2 %in% W))| (!(m1 %in% W))|(!(m2 %in% W))|(!(p2 %in% W))|(!(p1 %in% W))){D[n,]$FLAG='BAD'}
}
#tped format
#3       rs1000002       0       183635768       C       C       T       C       C       C
#3       rs1000003       0       98342907        A       A       A       A       A       A
#4       rs10000030      0       103374154       G       A       G       A       G       G
str=paste(D[n,]$CHR,D[n,]$uploaded_variant, "0",D[n,]$POS, D[n,]$C1, D[n,]$C2,D[n,]$M1,
          D[n,]$M2, D[n,]$P1, D[n,]$P2, sep = "  ");
print(str)

for(n in 1:N){
  c1=as.character(D[n,]$C1)
  c2=as.character(D[n,]$C2)
  m1=as.character(D[n,]$M1)
  m2=as.character(D[n,]$M2)
  ref=as.character(D[n,]$REF)
  alt=as.character(D[n,]$ALT)
  p1=as.character(D[n,]$P1)
  p2=as.character(D[n,]$P2)
  
  W=c(ref, alt )
##  if(!(c1 %in% W)){c1=complement(c1);D[n,]$C1=c1; print('complement')}
  if(!(c2 %in% W)){c2=complement(c2);D[n,]$C2=c2; print('complement')}
  if(!(m1 %in% W)){m1=complement(m1);D[n,]$M1=m1; print('complement')}
  if(!(m2 %in% W)){m2=complement(m2);D[n,]$M2=m2; print('complement')}
  if(!(p1 %in% W)){p1=complement(p1);D[n,]$P1=p1; print('complement')}
  if((p2 %in% W)&(p1 %in% W)&
     (m2 %in% W)&(m1 %in% W)&
     (c2 %in% W)&(c1 %in% W)){
  str=paste(D[n,]$CHR,D[n,]$uploaded_variant, "0",D[n,]$POS, D[n,]$C1, D[n,]$C2,D[n,]$M1,
            D[n,]$M2, D[n,]$P1, D[n,]$P2, sep = "  ");
  write.table(str,'trio.2.tped', quote=F, row.names = F, col.names = F, append = T);
}}

