#coral_reproduction_functions.R


get_protein_names = function(df, prot.col){
	ll=load("~/gitreps/coral_reproductive_evolution/datasets/ProteinTable.Rdata")
	ptable[,prot.col] = ptable$genbank.prot
	x=merge(df, ptable, by = prot.col, all.x = T)
	keep = append(colnames(df), 'protein.name')
	res=x[,keep]
	return(res)
	
}


make_sig_matrix = function(df, ptype, alpha){
	letters = c('g', 'm', 'p', 'q', 'v') #galaxia, montipora, porites, pocilloporids, all vertical
	pcols = paste(letters, bs.ptype, sep="_")
	pvals = df[,pcols]                               #subset the branch sites p-values
	sig = apply(pvals, 2, function(x) x<alpha)  
	return(sig)
}



test_cooccurance = function(sig.matrix, col1, col2){
	a0=sig.matrix[,col1]
	b0=sig.matrix[,col2]
	a=a0[!is.na(a0) & !is.na(b0)]
	b=b0[!is.na(a0) & !is.na(b0)]
	c = data.frame(a,b)
	sums = apply(c, 1, sum)
	tsums = table(sums)
	both = tsums['2']
	niether = tsums['0']
	single=c[sums==1,]
	aonly = sum(single[,'a'])
	bonly = sum(single[,'b'])
	r1 = c(both, aonly)
	r2 = c(bonly, niether)
	t = rbind(r1, r2)
	print(t)
	fisher.test(t, alternative='greater')
}





bs_sig_overlap = function(df, ptype){
	letters = c('g', 'm', 'p', 'q') #galaxia, montipora, porites, pocilloporids, all vertical
	pcols = paste(letters, bs.ptype, sep="_")
	pvals = df[,pcols]                               #subset the branch sites p-values
	sig = apply(pvals, 2, function(x) x<bs.alpha)      #set significance for branch sites tests
	tsig = apply(sig, 1, function(x) sum(x, na.rm=T))  #total significant
	table(tsig)                                        #numbers of significant lineages across orthologs
	bs.genes = rownames(pvals[tsig>=bs.min.lineage,])  #genes significant in enough lineages
	bso = length(bs.genes)
	print(paste(bso, "total genes fit these branch sites criteria"))
	return(bs.genes)
}







