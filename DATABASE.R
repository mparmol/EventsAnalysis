library("taxonomizr")
library("ggplot2")
library("ggtree")
library("ape")
library("stringr")
library("gtools")
library("phangorn")
library(ggnewscale)
library("cowplot")
library(data.table)
library(genbankr)
					
isEmpty <- function(x) { #This function checks if a data frame is empty or not
   return(length(x)==0)
}

paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

##########Preparation ONLY RUN ONE TIME WHEN YOU GET NEW DATA

### 1 - Taxonomy database

prepareDatabase('/storage/parras/Taxa_ncbi_NEW/accessionTaxa.sql') ##Modify this path manually. Remember to change it in the other script too

print("Taxonomy database ready")

### 2 - Get an original dataframe from all sequence taxonomy

file_aa_list=system('find /storage/shared/ncbi_bacteria_assembly/GCA/ -name "predicted-orfs-amino.fasta"', intern = TRUE)

for(o in seq(1,length(file_aa_list),500))
{
	system(paste("cat ",paste(file_aa_list[o:(o+499)],collapse=" ")," >> Resistance_genes_list.fa"))
}

tablita<-read.table("Resistance_genes_list.fa",stringsAsFactors = FALSE)
tablita2<-data.frame(tablita[ seq(1, nrow(tablita), by = 2),])
jolin<-as.data.frame(t(data.frame(strsplit(as.character(tablita2[,1]),"_"))))
jolin2<-as.character(jolin$V2) ####Para bacterias

id<-accessionToTaxa(jolin2,"/storage/parras/Taxa_ncbi/accessionTaxa.sql") ##Modify this path manually. Remember to change it in the other script too
tab_resul<-data.frame(getTaxonomy(id,"/storage/parras/Taxa_ncbi/accessionTaxa.sql")) ##Modify this path manually. Remember to change it in the other script too
tab_resul$species<-str_replace(tab_resul$species,"\\[","")
tab_resul$species<-str_replace(tab_resul$species,"\\]","")
tab_resul$species<-str_replace(tab_resul$species,"\\:","")
tab_resul$species<-str_replace(tab_resul$species,"\\(","")
tab_resul$species<-str_replace(tab_resul$species,"\\)","")
tab_resul$species<-str_replace_all(tab_resul$species,"\\'","")
write.table(tab_resul,"/home/parras/RESULTADOS_FINALES_NO_BORRAR/complete_taxa_file.txt",sep="\t") ##Modify this path manually. Remember to change it in the other script too

tablita[ seq(1, nrow(tablita), by = 2),]<-paste(">",jolin$V2,"_",jolin$V3,"-",gsub(" ", "_", as.character(tab_resul$species), fixed = TRUE),sep = '')
name_tablita<-data.frame(paste(jolin$V2,"_",jolin$V3,"-",gsub(" ", "_", as.character(tab_resul$species), fixed = TRUE),sep = ''))

write.table(sapply(strsplit(unique(as.character(name_tablita[,1])),"_"),"[[",1),"/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_acc_list.txt",row.names=F,col.names=F,quote=F,sep="\t") ##Modify this path manually. Remember to change it in the other script too

system("rm -rf Resistance_genes_list.fa")

print("Taxonomy tables ready")

### 3 - Get accession number and names for each analyzed sequence

system("find /storage/shared/ncbi_bacteria_assembly/GCA/ -name '*.fna' -exec grep '>' {} \\; > /home/parras/RESULTADOS_FINALES_NO_BORRAR/total_ids.txt") ##Modify this path manually. Remember to change it in the other script too

ids_total=data.frame(fread("/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_ids.txt",header = F,sep = "")) ##Modify this path manually. Remember to change it in the other script too
ids_total=cbind(sapply(strsplit(sapply(strsplit(x = ids_total[,1]," "),"[[",1),">"),"[[",2),ids_total)

tab_ids_totales=read.delim("/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_acc_list.txt",header=F) ##Modify this path manually. Remember to change it in the other script too

ids_total_saved=ids_total[ids_total[,1] %in% tab_ids_totales[,1],]

write.table(ids_total_saved,"/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_ids.txt",row.names=F,col.names=F,quote=F,sep="\t") ##Modify this path manually. Remember to change it in the other script too

system("rm -rf /home/parras/RESULTADOS_FINALES_NO_BORRAR/total_acc_list.txt") ##Modify this path manually. Remember to change it in the other script too

print("Get all ID genome list ready")

### 4 - Download all information related to each sequence, like isolation_source
##Get bact info

system("wget ftp.ncbi.nlm.nih.gov/genbank/")

index_gb=read.delim("index.html")
file_list=sapply(strsplit(sapply(strsplit(as.character(index_gb[grep("gbbct",index_gb[,1]),1]),">"),"[[",1),"="),"[[",2)

for(o in 1:length(file_list))
{
	system(paste("wget ftp.ncbi.nlm.nih.gov/genbank/",file_list[o],sep=""))
}

system("gzip -d *gz")

files_totales=list.files("./")

res_isol=NULL
e=1

for(i in 1:length(files_totales))
{
  if(!isEmpty(grep(".seq",files_totales[i])))
  {
    p=10
    system(paste('sed -e "1,',p,'d" < ',files_totales[i],' > file_filtered.gb',sep=""))

    while(isEmpty(system('head -1 file_filtered.gb | grep "LOCUS"',intern=T)))
    {
      p=p-1
      system(paste('sed -e "1,',p,'d" < ',files_totales[i],' > file_filtered.gb',sep=""))   
    }
		
    system('grep "/chromosome" -v file_filtered.gb > file_filtered2.gb')
    
    system("awk -v n=1 '/^\\/\\//{close(\"out_\"n);n++;next} {print > \"out_\"n}' file_filtered2.gb")
    
    files_totales2=list.files("./")
    
    pg.list <- lapply(files_totales2[grep("out_",files_totales2)], function(x){parseGenBank(x,ret.seq = FALSE)})
    
    aux_isol=data.frame(matrix(nrow = length(pg.list),ncol=11))
    
    for(w in 1:length(pg.list))
    {
      if(!isEmpty(paste(pg.list[[w]]$LOCUS,collapse=""))){aux_isol[w,1]=paste(pg.list[[w]]$LOCUS,collapse="")}
      if(!isEmpty((pg.list[[w]]$VERSION)[[1]])){aux_isol[w,2]=(pg.list[[w]]$VERSION)[[1]]}
      if(!isEmpty(pg.list[[w]]$FEATURES$`1`$isolation_source[1])){aux_isol[w,5]=pg.list[[w]]$FEATURES$`1`$isolation_source[1]}
      if(!isEmpty(pg.list[[w]]$FEATURES$`1`$country[1])){aux_isol[w,6]=pg.list[[w]]$FEATURES$`1`$country[1]}
      if(!isEmpty(pg.list[[w]]$FEATURES$`1`$collection_date[1])){aux_isol[w,10]=pg.list[[w]]$FEATURES$`1`$collection_date[1]}
      if(!isEmpty(pg.list[[w]]$FEATURES$`1`$host[1])){aux_isol[w,11]=pg.list[[w]]$FEATURES$`1`$host[1]}
    }
    
    res_isol=rbind(res_isol,aux_isol)
    
    system("rm -rf out* file_filtered.gb file_filtered2.gb")
  }
  
  write.table(res_isol,"/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt",col.names=F,row.names=F,sep="\t",quote=F) ##Modify this path manually. Remember to change it in the other script too

}

write.table(res_isol,"/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt",col.names=F,row.names=F,sep="\t",quote=F) ##Modify this path manually. Remember to change it in the other script too

system("rm -rf *seq index.html")

##Get WGS info

e=1

system("wget https://ftp.ncbi.nlm.nih.gov/genbank/wgs/")

system('sed -e "1,8d" < index.html > index2.html')

system("cut -d'/' -f 1 index2.html | grep -v '\\.' | cut -d'\"' -f 2 | sort > index3.html")

tab_wgs=read.delim("index3.html")

for(i in 1:length(tab_wgs[,1]))
{
	if(nchar(as.character(tab_wgs[i,1]))==3)
	{
		system(paste("wget https://ftp.ncbi.nlm.nih.gov/genbank/wgs/",tab_wgs[i,1],"/",sep=""))
		system("grep 'gbff' index.html.1 | grep -v 'mstr' | cut -d'>' -f 2 | cut -d'<' -f 1 > index_sub.html")
		
		tab_wgs_2=read.delim("index_sub.html")

		for(o in 1:length(tab_wgs_2[,1]))
		{
			system(paste("wget https://ftp.ncbi.nlm.nih.gov/genbank/wgs/",tab_wgs[i,1],"/",tab_wgs_2[o,1],sep=""))
		}

		system("gzip -d *gz")

		files_totales=list.files("./")
   

		system("rm -rf index.html.1")
 
		for(k in 1:length(files_totales))
		{
  			if(!isEmpty(grep("gbff",files_totales[k])))
  			{
		   
    				#system(paste('sed -e "1,10d" < ',files_totales[k],' > file_filtered.gb',sep=""))
	 
  
    				system(paste('grep "/chromosome" -v ',files_totales[k],' > file_filtered2.gb'))
    
    				system("awk -v n=1 '/^\\/\\//{close(\"out_\"n);n++;next} {print > \"out_\"n}' file_filtered2.gb")
    
    				files_totales2=list.files("./")
    
    				pg.list <- lapply(files_totales2[grep("out_",files_totales2)], function(x){parseGenBank(x,ret.seq = FALSE)})
    
    				aux_isol=data.frame(matrix(nrow = length(pg.list),ncol=11))
    
    				for(w in 1:length(pg.list))
    				{
      					if(!isEmpty(paste(pg.list[[w]]$LOCUS,collapse=""))){aux_isol[w,1]=paste(pg.list[[w]]$LOCUS,collapse="")}
      					if(!isEmpty((pg.list[[w]]$VERSION)[[1]])){aux_isol[w,2]=(pg.list[[w]]$VERSION)[[1]]}
      					if(!isEmpty(pg.list[[w]]$FEATURES$`1`$isolation_source[1])){aux_isol[w,5]=pg.list[[w]]$FEATURES$`1`$isolation_source[1]}
      					if(!isEmpty(pg.list[[w]]$FEATURES$`1`$country[1])){aux_isol[w,6]=pg.list[[w]]$FEATURES$`1`$country[1]}
      					if(!isEmpty(pg.list[[w]]$FEATURES$`1`$collection_date[1])){aux_isol[w,10]=pg.list[[w]]$FEATURES$`1`$collection_date[1]}
      					if(!isEmpty(pg.list[[w]]$FEATURES$`1`$host[1])){aux_isol[w,11]=pg.list[[w]]$FEATURES$`1`$host[1]}
    				}
    
    				res_isol=rbind(res_isol,aux_isol)
    
    				system("rm -rf out* file_filtered2.gb")
  			}
  
  			write.table(res_isol,"/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt",col.names=F,row.names=F,sep="\t",quote=F) ##Modify this path manually. Remember to change it in the other script too
		}
 
		write.table(res_isol,"/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt",col.names=F,row.names=F,sep="\t",quote=F) ##Modify this path manually. Remember to change it in the other script too
		system("rm -rf *gbff")
		
	}
}

print("Isolation source database ready")

### 5 - Download ResFinder database file (nucleotide format db).

system("git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git")
system("cd resfinder_db")
system("cat *fsa > complete_db")
system("makeblastdb -in complete_db -dbtype nucl -out /storage/parras/ResFinder_db/AntibioticDatabase") ##Modify this path manually. Remember to change it in the other script too
print("ResFinder database file ready")
