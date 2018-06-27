#get_gene_ids.R
#This script gathers gene descriptions from the swissprot annotations

#UPLOAD ANNOTATIONS
adat = read.table("~/gitreps/convergent_evo_coral/ortholog_tables/singleCopyAnnotations.tsv", header = T)[,c(1,2)]
colnames(adat) = c('seq', 'swissProt_hits')

#FORMAT ANNOATIONS
swissProt = c()
ortho = c()
name = c()
for (i in 1:nrow(adat)){
	row=adat[i,]
	annots = unlist(strsplit(as.character(row$swissProt_hits), ';'))
	sp = unlist(sapply(annots, function(x) strsplit(x, "|", fixed=T)[[1]][1]))
	id = unlist(sapply(annots, function(x) strsplit(x, "|", fixed=T)[[1]][2]))
	seqs = rep(as.character(row$seq), length(sp))
	ortho = append(ortho, seqs)
	swissProt = append(swissProt, sp)
	name = append(name, id)
}
adf = data.frame(ortho, swissProt, name)
db = sapply(as.character(adf$name), function(x) strsplit(x, "_")[[1]][2])
egn = sapply(as.character(adf$name), function(x) strsplit(x, "_")[[1]][1])
adf$egn = egn
adf$db=db
head(adf)
length(unique(adf$db))


#----- GET GENE DESCRIPTIONS -----#
#unfortunately, the descriptions have to be downloaded by species
#and these are split among many different ones
#will add them iteratively from each species

#ADD HUMAN ANNOTATIONS
library("biomaRt")

add.species = function(dataset.name, abbreviation){
	embl = useMart("ensembl", dataset=dataset.name) #get the human dataset
	sub = adf[adf$db== abbreviation & !is.na(adf$db),]
	spSet = unique(sub$swissProt)
	gnames.sp = getBM(attributes = c('uniprotswissprot', 'external_gene_name','description'), filters = c('uniprotswissprot'), values = spSet, mart = embl)
	
	get.description = function(x){
	g=gnames.sp[gnames.sp$uniprotswissprot==x, 'description']
	if (length(g) < 1){
		return('none')
	}
	else{
		return(g[1])
		}
	}
	
	des =unlist(sapply(sub$swissProt, function(x) get.description(x)))
	res = data.frame(sub$swissProt, des)
	colnames(res) = c('swissProt2', 'description')
	sub2 = cbind(sub, res)
	print("Everything Match?")
	print(sum(sub2$swissProt== sub2$swissProt2)==nrow(sub2))
	sub2$swissProt2<-NULL
	return(sub2)
}

combine.annotations = function(annotated, new){
	good = annotated[annotated$description!="none",]
	to.add = new[!new$ortho %in% good$ortho & new$description != "none",]
	combined = rbind(good, to.add)
	td = length(unique(combined$ortho))
	ta = length(unique(adf$ortho))
	pct = round(td / ta, digits=4) * 100
	print("Percentage of orthologs annotated:")
	print(paste(pct, "%", sep = ''))
	combined=combined[!duplicated(combined$ortho),]
	return(combined)
}

#look at available datasets
mart = useMart("ensembl")
mart.dats = listDatasets(mart)
head(mart.dats)
mart.dats$dataset[grep("ele", mart.dats$dataset)]

#start annotations off with human
x = add.species("hsapiens_gene_ensembl", "HUMAN")
annotated = x
td = length(unique(annotated$ortho))
ta = length(unique(adat$seq))
pct = round(td / ta, digits=4) * 100
paste(pct, "%", sep = '')

#now add remaining species
x = add.species("mmusculus_gene_ensembl", "MOUSE") #human
annotated = combine.annotations(annotated, x)
x = add.species("rnorvegicus_gene_ensembl", "RAT") #rat
annotated = combine.annotations(annotated, x)
x = add.species("btaurus_gene_ensembl", "BOVIN") #cow
annotated = combine.annotations(annotated, x)
x = add.species("drerio_gene_ensembl", "DANRE") #zebrafish
annotated = combine.annotations(annotated, x)
x = add.species("ggallus_gene_ensembl", "CHICK") #chicken
annotated = combine.annotations(annotated, x)
x = add.species("dmelanogaster_gene_ensembl", "DROME") #fruit fly
annotated = combine.annotations(annotated, x)
x = add.species("sscrofa_gene_ensembl", "PIG") #pig
annotated = combine.annotations(annotated, x)


#assemble descriptions
head(annotated)
head(adf)
no.d = adf[!adf$ortho %in% annotated$ortho,]
head(no.d)
adf$description = 'none'
final = rbind(annotated, adf)
final = final[!duplicated(final$ortho),]
head(final)
write.table(final, file="~/gitreps/convergent_evo_coral/ortholog_tables/singleCopyAnnotationsWithDescriptions.tsv")


