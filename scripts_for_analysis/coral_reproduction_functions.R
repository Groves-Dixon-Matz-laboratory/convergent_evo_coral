#coral_reproduction_functions.R


readin_branch_sites_lrt_results = function(file_name, clade){
	dat = read.table(file_name, sep= "\t", header = T, quote="", stringsAsFactors=F)
	dat$clade = clade
	return(dat)
}

#--- filter_bs()
# #returns the list of genes that pass the branch site test filters:
# bs.ptype            #type of p value to use (adj.p or pvalues)
# bs.alpha            #significance cutoff for chosen p value
# bs.min.lineage = 1  #number of target lineages that must be significant to keep gene

filter_bs = function(dat){
	sig <- dat %>% filter(dat[,bs.ptype] < bs.alpha)
	counts = table(sig$ortholog)
	keep = names(counts[counts>=bs.min.lineage])
	filt = dat[dat$ortholog %in% keep,]
	nGenes0 = length(unique(dat$ortholog))
	nGenesF = length(unique(filt$ortholog))
	propFiltered = nGenesF / nGenes0
	pctFiltered = round(propFiltered, digits=3)*100
	print("Filtering...")
	print("Clades:")
	print(unique(sig$clade))
	print(paste("Using p-value type:", bs.ptype))
	print(paste("with cutoff <", bs.alpha))
	print(paste("For a gene to pass, number of significant lineages must be at least", bs.min.lineage))
	print(paste(nGenes0, 'genes before filtering'))
	print(paste(nGenesF, 'genes after filtering'))
	print(paste(pctFiltered, '% genes passing', sep=''))
	return(filt$ortholog)
}


filter_flags = function(dat){
	filt0 <- dat %>% 
		filter(posterior >= bs.site.post.cutoff) %>%
		filter(bs_passable)
	counts = table(filt0$site)
	keep = names(counts[counts>= min.flagged.tests])
	filt = filt0[filt0$site %in% keep,]
	nSubs0 = nrow(dat)
	nSubsF = nrow(filt)
	propSubsFiltered = nSubsF / nSubs0
	pctSubsFiltered = round(propSubsFiltered, digits=3)*100
	nGenes0 = length(unique(dat$ortholog))
	nGenesF = length(unique(filt$ortholog))
	propFiltered = nGenesF / nGenes0
	pctFiltered = round(propFiltered, digits=3)*100
	print("Filtering...")
	print("Clades:")
	print(unique(filt$clade))
	print(paste("substitutions must have posterior >=", bs.site.post.cutoff))
	print(paste(nSubs0, 'substitutions before filtering'))
	print(paste(nSubsF, 'substitutions after filtering'))
	print(paste(pctSubsFiltered, '% substitutions passing', sep=''))
	print('.')
	print(paste("gene must have at least one site with posterior >=", bs.site.post.cutoff))
	print(paste(nGenes0, 'genes before filtering'))
	print(paste(nGenesF, 'genes after filtering'))
	print(paste(pctFiltered, '% genes passing', sep=''))
	return(filt$site)
}


filter_convergence = function(dat0){
	#get convergent only
	dat <- dat0 %>% filter(genConvergent)
	
	#apply posterior filtering for convergence
	dat$site = paste(dat$ortholog, dat$b1pos, sep='_')
	sub <- dat %>% 
	filter(b1post >= anc.post.cutoff) %>%
	filter(b2post >= anc.post.cutoff)

	#subset for sites with at least min.convergent.changes
	tdf <- sub %>% group_by(site) %>%
	mutate(nClades = length(unique(append(as.character(b1Clade), as.character(b2Clade))))) %>%
	select(c(ortholog, b1pos, site, nClades))
	keep = tdf$site[tdf$nClades >= min.convergent.changes]
	filt <- sub %>% filter(site %in% keep)
	
	#count filtering subs
	nSubs0 = length(unique(dat$site))
	nSubs1 = length(unique(sub$site)) #n unique subs after posterior filtering
	propSubsFiltered = nSubs1 / nSubs0
	pctSubsFiltered1 = round(propSubsFiltered, digits=3)*100
	
	#count filtering genes
	nGenes0 = length(unique(dat$ortholog))
	nGenes1 = length(unique(sub$ortholog))
	propGenesFiltered = nGenes1 / nGenes0
	pctGenesFiltered1 = round(propGenesFiltered, digits=3)*100
	
	#count filtering subs 2
	nSubs2 = length(unique(filt$site))
	propSubsFiltered = nSubs2 / nSubs0
	pctSubsFiltered2 = round(propSubsFiltered, digits=3)*100
	
	#count filtering genes 2
	nGenes2 = length(unique(filt$ortholog))
	propGenesFiltered = nGenes2 / nGenes0
	pctGenesFiltered2 = round(propGenesFiltered, digits=3)*100
	
	#report filtering results
	print("Filtering on posteriors...")
	print("Clades:")
	print(unique(dat$clade))
	print(paste("convergence calls must have posterior >=", anc.post.cutoff))
	print(paste(nSubs0, 'substitutions before filtering'))
	print(paste(nSubs1, 'substitutions after filtering'))
	print(paste(pctSubsFiltered1, '% substitutions passing', sep=''))
	print('.')
	print(paste("genes with at least one convergence with posterior >=", anc.post.cutoff))
	print(paste(nGenes0, 'genes before filtering'))
	print(paste(nGenes1, 'genes after filtering'))
	print(paste(pctGenesFiltered1, '% genes passing', sep=''))
	print('...')
	print('...')
	print("Filtering on number of convergent changes...")
	print(paste("site must have N convergent changes >=", min.convergent.changes))
	print(paste(nSubs0, 'substitutions before filtering'))
	print(paste(nSubs2, 'substitutions after filtering'))
	print(paste(pctSubsFiltered2, '% substitutions passing', sep=''))
	print('.')
	print(paste("genes with at least one site with N convergence changes >=", anc.post.cutoff))
	print(paste(nGenes0, 'genes before filtering'))
	print(paste(nGenes2, 'genes after filtering'))
	print(paste(pctGenesFiltered2, '% genes passing', sep=''))
	return(filt$site)
}

#Functions to build column that states:

# #Get set of unique bs filter-passing genes for each clade
# gather_pos_genes_for_flag_pass = function(x){
	# csub <- bs_filt %>%
		# filter(clade==x & p.values < bs.gene.pvalue.cutoff)
	# return(unique(csub$ortholog))
# }

# cross_check_bs_passable = function(x){
	# c = x['clade']
	# gene=x['ortholog']
	# geneSet = passableGeneList[[c]]
	# return(gene %in% geneSet)
# }




#Get set of unique bs filter-passing genes for each clade
gather_pos_genes = function(x){
	csub <- bs_filt %>%
		filter(clade==x & bs_pass)
	return(unique(csub$ortholog))
}

#Crosscheck those against genes convergence events occured in
#for a given overlapping substitition, 
#is the gene the substitition was found in
#among the genes with significant branch site results
#for clade1 of the convergence event, or for the 'All'
#branch site test that matches the phenotype for clade1

cross_check_bs_clade1 = function(x){
	c = x['b1Clade']
	gene=x['ortholog']
	if (grepl('Anti', c)){
		all = posGeneList[['AntiAllTarget']]
	} else {
		all = posGeneList[['AllTarget']]
	}
	geneSet = append(all, posGeneList[[c]])
	return(gene %in% geneSet)
}

#same as above but for clade 2
cross_check_bs_clade2 = function(x){
	c = x['b2Clade']
	gene=x['ortholog']
	if (grepl('Anti', c)){
		all = posGeneList[['AntiAllTarget']]
	} else {
		all = posGeneList[['AllTarget']]
	}
	geneSet = append(all, posGeneList[[c]])
	return(gene %in% geneSet)
}



#Function to build column that states:
#for a given overlapping substitition, 
#is the site among those flagged for positive
#selection using BEB in the branch-site test for
#clade1, or in the 'All' branch-site test for 
#the phenotype matching clade1
cross_check_flags_clade1 = function(x){
	c=x['b1Clade']
	site=x['site']
	if (grepl('Anti', c)){
		all = flagSiteList[['AntiAllTarget']]
	} else {
		all = flagSiteList[['AllTarget']]
	}
	siteSet = append(all, flagSiteList[[c]])
	return(site %in% siteSet)
}

#same for clade2
cross_check_flags_clade2 = function(x){
	c=x['b2Clade']
	site=x['site']
	if (grepl('Anti', c)){
		all = flagSiteList[['AntiAllTarget']]
	} else {
		all = flagSiteList[['AllTarget']]
	}
	siteSet = append(all, flagSiteList[[c]])
	return(site %in% siteSet)
}


#check if a site passes all conditions 
cross_check_convflags = function(x){
	cpair=x['cp']
	gene = x['ortholog']
	geneSet = flagConvGeneList[[cpair]]
	return(gene %in% geneSet)
}


#check if a site passes all conditions 
cross_check_all = function(x){
	cpair=x['cp']
	gene = x['ortholog']
	geneSet = allEvList[[cpair]]
	return(gene %in% geneSet)
}


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
	abrvs = c('g', 'm', 'p', 'q', 'v') #galaxia, montipora, porites, pocilloporids, all vertical
	pcols = paste(abrvs, bs.ptype, sep="_")
	pvals = df[,pcols]                               #subset the branch sites p-values
	sig = apply(pvals, 2, function(x) x<bs.alpha)      #set significance for branch sites tests
	tsig = apply(sig, 1, function(x) sum(x, na.rm=T))  #total significant
	table(tsig)                                        #numbers of significant lineages across orthologs
	bs.genes = rownames(pvals[tsig>=bs.min.lineage,])  #genes significant in enough lineages
	bso = length(bs.genes)
	print(paste(bso, "total genes fit these branch sites criteria"))
	return(bs.genes)
}







