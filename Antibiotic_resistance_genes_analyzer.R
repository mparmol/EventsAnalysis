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
library("adephylo")
library("randomcoloR")
library("stringr")

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

############ Main

#model_list=c("mph.hmm","tet_enzyme.hmm","aac3_class1.hmm","aac6p_class3.hmm","aph2b.hmm","aac2p.hmm","class_B_1_2.hmm","class_D_1.hmm","aac6p_class1.hmm","erm_typeF.hmm","qnr.hmm","class_C.hmm","aac3_class2.hmm","tet_rpg.hmm","aph6.hmm","class_D_2.hmm","tet_efflux.hmm","aph3p.hmm","class_A.hmm","class_B_3.hmm","aac6p_class2.hmm","B1.hmm")
model_list=c("erm_typeF.hmm")

file_aa_list=system('find /storage/shared/ncbi_bacteria_assembly/GCA/ -name "predicted-orfs-amino.fasta"', intern = TRUE)

print("Files read")

total_ids=data.frame(fread("/home/parras/RESULTADOS_FINALES_NO_BORRAR/total_ids.txt",header = F,sep = "\t"))


print("Total IDs read")

for(j in 1:length(model_list))
{

system(paste('mkdir ',model_list[j]))

sublist=file_aa_list[grep(model_list[j],file_aa_list)]

for(o in seq(1,length(sublist),500))
{
	system(paste("cat ",paste(sublist[o:(o+499)],collapse=" ")," >> Resistance_genes_list.fa"))
}

#system(paste('cat /storage/shared/ncbi_bacteria_assembly/GCA/*/*/*/*/',model_list[j],'/predictedGenes/predicted-orfs-amino.fasta > Resistance_genes_list.fa',sep="")) ##Get all genes from a certain class

sequence_info<-read.table("Resistance_genes_list.fa",stringsAsFactors = FALSE)
sequence_info2<-data.frame(sequence_info[ seq(1, nrow(sequence_info), by = 2),])
sequence_info_j<-as.data.frame(t(data.frame(strsplit(as.character(sequence_info2[,1]),"_"))))
sequence_info_j2<-as.character(sequence_info_j$V2) ####Para bacterias

id<-accessionToTaxa(sequence_info_j2,"/storage/parras/Taxa_ncbi/accessionTaxa.sql")
tab_resul<-data.frame(getTaxonomy(id,"/storage/parras/Taxa_ncbi/accessionTaxa.sql"))
tab_resul$species<-str_replace(tab_resul$species,"\\[","")
tab_resul$species<-str_replace(tab_resul$species,"\\]","")
tab_resul$species<-str_replace(tab_resul$species,"\\:","")
tab_resul$species<-str_replace(tab_resul$species,"\\(","")
tab_resul$species<-str_replace(tab_resul$species,"\\)","")
tab_resul$species<-str_replace_all(tab_resul$species,"\\'","")
#write.table(sequence_info_j,"sequence_info_j.txt",sep="\t")

sequence_info[ seq(1, nrow(sequence_info), by = 2),]<-paste(">",sequence_info_j$V2,"_",sequence_info_j$V3,"-",gsub(" ", "_", as.character(tab_resul$species), fixed = TRUE),sep = '')
name_sequence_info<-data.frame(paste(sequence_info_j$V2,"_",sequence_info_j$V3,"-",gsub(" ", "_", as.character(tab_resul$species), fixed = TRUE),sep = ''))
row.names(name_sequence_info)<-name_sequence_info[,1]
row.names(tab_resul)<-row.names(name_sequence_info)

#write.table(tab_resul,"prueba_a_ver_ids_mod.txt",sep="\t",quote = FALSE)

write.table(sequence_info,"Resistance_genes_list_species.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)


#####Complet_taxa

taxonomy_db=data.frame(fread("/home/parras/RESULTADOS_FINALES_NO_BORRAR/complete_taxa_file.txt",header = F))

id_mods_=tab_resul

id_mods_[,7]=as.character(id_mods_[,7])

for(o in 1:length(id_mods_[,1]))
{

	if(is.na(id_mods_[o,2]))
	{
		
		taxonomic_info=strsplit(as.character(total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],])[2]," ")[[1]][2]

		if(!is.na(taxonomic_info))
		{

		if(!isEmpty(grep("TPA_asm",taxonomic_info)))
		{
			#print(taxonomic_info)
			taxonomic_info=strsplit(as.character(total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],])[2]," ")[[1]][3]
			#print(taxonomic_info)
		}

		if(!isEmpty(grep(taxonomic_info,taxonomy_db[,7])))
		{
			id_mods_[o,]=as.vector(c(na.omit(taxonomy_db[taxonomy_db[,7]==taxonomic_info,])[1,2:7],gsub(">","",total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],][2])))
		}else if(!isEmpty(grep(taxonomic_info,taxonomy_db[,6])))
		{
			id_mods_[o,]=as.vector(c(na.omit(taxonomy_db[taxonomy_db[,6]==taxonomic_info,])[1,2:6],"NA",gsub(">","",total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],][2])))

		}else if(!isEmpty(grep(taxonomic_info,taxonomy_db[,5])))
		{
			id_mods_[o,]=as.vector(c(na.omit(taxonomy_db[taxonomy_db[,5]==taxonomic_info,])[1,2:5],"NA","NA",gsub(">","",total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],][2])))

		}else if(!isEmpty(grep(taxonomic_info,taxonomy_db[,4])))
		{
			id_mods_[o,]=as.vector(c(na.omit(taxonomy_db[taxonomy_db[,4]==taxonomic_info,])[1,2:4],"NA","NA","NA",gsub(">","",total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],][2])))

		}else if(!isEmpty(grep(taxonomic_info,taxonomy_db[,3])))
		{
			id_mods_[o,]=as.vector(c(na.omit(taxonomy_db[taxonomy_db[,3]==taxonomic_info,])[1,2:3],"NA","NA","NA","NA",gsub(">","",total_ids[total_ids[,1]==strsplit(rownames(id_mods_)[o],"_")[[1]][1],][2])))
	
		}

		}


	
	}


	
	

}

id_mods_[grep("Candidatus",id_mods_[,2]),2:7]=NA
id_mods_[grep("candidate",id_mods_[,2]),2:7]=NA #Candidatus entries are removed from the taxonomy infor, they won't create events

write.table(id_mods_,"Taxonomy_info_dataset.txt",sep="\t",quote=F)


	system("clustalo -i Resistance_genes_list_species.fasta -o Alignment_file")
	system("~/.conda/envs/marcosEnv/bin/fasttree Alignment_file > Tree_file")
	system("~/.conda/envs/marcosEnv/bin/tblastn -max_hsps 1 -max_target_seqs 1 -query Resistance_genes_list_species.fasta -outfmt 6 -db /storage/parras/ResFinder_db/AntibioticDatabase -out Blast_Antibiotic_analysis")

tab_resul<-id_mods_ #Here we read the taxonomic information table we are going to use to evaluate events

name_sequence_info<-data.frame(rownames(tab_resul))

name_sequence_info$phylum<-tab_resul$phylum

tree_l<-read.tree("Tree_file") #The unrooted tree is loaded this way
tree_l<-makeNodeLabel(tree_l) #With this funtion we give name to each of the tree nodes
write.tree(tree_l,"Tree_file")

row.names(name_sequence_info)<-rownames(tab_resul) #We create an auxiliary table to save all the information we are processing during the analysis, with taxonomic information, blast best hit...
name_sequence_info2<-data.frame(name_sequence_info[,-1])
row.names(name_sequence_info2)<-row.names(name_sequence_info)

tree<-tree_l

####With this code block we are able to get all the information from the nodes and leaves, their ID and given name

node_labels_in_edge <- tree$node.label[tree$edge[,1]-Ntip(tree)]
tips_nodes <- tree$edge[,2]

select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}

edge_table <- data.frame(
  "parent" = tree$edge[,1],
  "par.name" = sapply(tree$edge[,1], select.tip.or.node, tree = tree),
  "child" = tree$edge[,2],
  "chi.name" = sapply(tree$edge[,2], select.tip.or.node, tree = tree)
)

edge_table$parent=as.numeric(edge_table$parent)
edge_table$child=as.numeric(edge_table$child)
edge_table$par.name=as.character(edge_table$par.name)
edge_table$chi.name=as.character(edge_table$chi.name)
####


name_sequence_info3<-name_sequence_info2 ####Auxiliary table only for phyla between node comparations

nodes_list=unique(edge_table$par.name) #We get all possible node names from the table
output_results=data.frame(matrix(ncol = 10))
output_results_pairwise=data.frame(matrix(ncol = 3))
e=1

for(i in 1:length(nodes_list)) #This loop is the core of the script. It will go over all nodes form last to first one, comparing the taxonomy of the sequeneces that belongs to each of the branches
{
  if(!sort(str_detect(edge_table[edge_table$child %in% Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]],4],"Node."),decreasing = T)[1]) ###This checks for a 100% similarity sequences nodes that could be composed of genomes belonging to different phyla 
  {
    child_node_levels<-levels(factor(as.character(name_sequence_info3[rownames(name_sequence_info3) %in% edge_table[edge_table$child %in% Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]],4],1])))
    
    if((length(child_node_levels)==2 && !('Unknown' %in% child_node_levels)) || length(child_node_levels)>2)
    {
      output_results[e,]=c(nodes_list[length(nodes_list)-i+1],paste(child_node_levels,collapse = ";"),NA,paste(edge_table[edge_table$child %in% Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]],4],collapse = ";"),NA,length(edge_table[edge_table$child %in% Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]],4]),NA,NA,NA,NA) #Within this line we save the node info of events
      output_results_pairwise[e,]<-c(nodes_list[length(nodes_list)-i+1],"NA","NA")
	e=e+1
    }
  }else
  {### In this block we check that a node contains only leaves and their phyla are different among them
    
    descendant_nodes=Descendants(tree,edge_table[edge_table$par.name==nodes_list[length(nodes_list)-i+1],1],type = "child")[[1]] #Knowing which are the children of each node we can get their taxonomy and compare it
    descendant_nodes1=Descendants(tree,descendant_nodes[1])[[1]]
    descendant_nodes2=Descendants(tree,descendant_nodes[2])[[1]]
    
    vector1=NULL #Here we save the leaves names from each descendant
    vector2=NULL
    
    for(o in 1:length(descendant_nodes1))
    {
      vector1=c(vector1,edge_table[edge_table$child==descendant_nodes1[o],4])  
    }
    
    for(o in 1:length(descendant_nodes2))
    {
      vector2=c(vector2,edge_table[edge_table$child==descendant_nodes2[o],4])  
    }
    
    vector1=sort(vector1) #Leaves names are sorted here
    vector2=sort(vector2)
    
    leaves_levels1<-levels(factor(as.character(name_sequence_info3[rownames(name_sequence_info3) %in% vector1,1])))
    leaves_levels2<-levels(factor(as.character(name_sequence_info3[rownames(name_sequence_info3) %in% vector2,1])))
    
    ################################
  
    if(length(descendant_nodes)==3) #This is only in case we are working with tripartite trees, so we can get the info for the three branches
    {
      descendant_nodes3=Descendants(tree,descendant_nodes[3])[[1]]
      
      vector3=NULL
      
      for(o in 1:length(descendant_nodes3))
      {
        vector3=c(vector3,edge_table[edge_table$child==descendant_nodes3[o],4])  
      }
      
      leaves_levels3<-levels(factor(as.character(name_sequence_info3[rownames(name_sequence_info3) %in% vector3,1])))
      
    }
	
    if(length(descendant_nodes)==3 && (leaves_levels1 != "Unknown" && leaves_levels2 != "Unknown" && leaves_levels3 != "Unknown" && length(unique(c(leaves_levels1,leaves_levels2,leaves_levels3)))>sort(c(length(leaves_levels1),length(leaves_levels2),length(leaves_levels3)),decreasing = T)[1]))
    {
		output_results[e,]=c(nodes_list[length(nodes_list)-i+1],paste(leaves_levels1,collapse = ";"),paste(leaves_levels2,collapse = ";"),paste(leaves_levels3,collapse = ";"),paste(vector1,collapse = ";"),paste(vector2,collapse = ";"),paste(vector3,collapse = ";"),length(vector1),length(vector2),length(vector3)) #Within this line we save the node info of events
        
		if(length(vector1)<=5000 & length(vector2)<=5000 & length(vector3)<=5000) #Here we calculate the inter-patristic distance within sequences that belongs to different branches, only for thoses that have 5000 or less leaves
        {
          sub_tree<-keep.tip(tree,c(vector1,vector2,vector3))
          
          pairwise_all<-data.frame(as.matrix(distTips(sub_tree,method = "patristic")))
          colnames(pairwise_all)<-rownames(pairwise_all)
          pairwise_all$leaf<-rownames(pairwise_all)
          pairwise_all_melt<-melt(pairwise_all,id.vars = "leaf")
          pairwise_all_melt$combination<-paste(pairwise_all_melt$leaf,pairwise_all_melt$variable,sep = "")
          
          
          leaves_combinations=data.frame(expand.grid(vector1,vector2,vector3)) #We perform the combinations and save them as characters
          leaves_combinations$Var1=as.character(leaves_combinations$Var1)
          leaves_combinations$Var2=as.character(leaves_combinations$Var2)
          leaves_combinations$Var3=as.character(leaves_combinations$Var3)
          values_j=NULL #Within this vector we will save all distance values from the pairwise of all leaves
          
          combinations_parse<-paste(leaves_combinations$Var1,leaves_combinations$Var2,leaves_combinations$Var3,sep = "")
          values_j<-pairwise_all_melt[pairwise_all_melt$combination %in% combinations_parse,3]
          
          output_results_pairwise[e,]<-c(nodes_list[length(nodes_list)-i+1],min(values_j),max(values_j))
        }else
        {
          output_results_pairwise[e,]<-c(nodes_list[length(nodes_list)-i+1],"NA","NA")
        }
        e=e+1  
      }else if(length(descendant_nodes)==2 && (leaves_levels1 != "Unknown" && leaves_levels2 != "Unknown" && length(unique(c(leaves_levels1,leaves_levels2)))>sort(c(length(leaves_levels1),length(leaves_levels2)),decreasing = T)[1]))
      {
  
       output_results[e,]=c(nodes_list[length(nodes_list)-i+1],paste(leaves_levels1,collapse = ";"),paste(leaves_levels2,collapse = ";"),paste(vector1,collapse = ";"),paste(vector2,collapse = ";"),length(vector1),length(vector2),NA,NA,NA) #Within this line we save the node info of events
          
       if(length(vector1)<=5000 & length(vector2)<=5000)
       {
        sub_tree<-keep.tip(tree,c(vector1,vector2))
            
        pairwise_all<-data.frame(as.matrix(distTips(sub_tree,method = "patristic")))
        colnames(pairwise_all)<-rownames(pairwise_all)
        pairwise_all$leaf<-rownames(pairwise_all)
        pairwise_all_melt<-melt(pairwise_all,id.vars = "leaf")
        pairwise_all_melt$combination<-paste(pairwise_all_melt$leaf,pairwise_all_melt$variable,sep = "")
            
            
        leaves_combinations=data.frame(expand.grid(vector1,vector2)) #We perform the combinations and save them as characters
        leaves_combinations$Var1=as.character(leaves_combinations$Var1)
        leaves_combinations$Var2=as.character(leaves_combinations$Var2)
        values_j=NULL 
            
        combinations_parse<-paste(leaves_combinations$Var1,leaves_combinations$Var2,sep = "")
        values_j<-pairwise_all_melt[pairwise_all_melt$combination %in% combinations_parse,3]
            
        output_results_pairwise[e,]<-c(nodes_list[length(nodes_list)-i+1],min(values_j),max(values_j))
       }else
       {
		output_results_pairwise[e,]<-c(nodes_list[length(nodes_list)-i+1],"NA","NA")
       }
       e=e+1  
      }  
  }
}

#write.table(output_results_pairwise,"output_results_pairwise.txt",sep = "\t",row.names = F,col.names = F,quote = F)

Nodes_to_color<-unique(edge_table[edge_table$par.name %in% output_results$X1,1])

output_results[,11]=NA
output_results[,12]=NA

################## Here we load the Blast result computed before

output_results_blast<-read.delim("Blast_Antibiotic_analysis",sep="\t",header=F)

name_sequence_info2[,2]<-output_results_blast[,2]
name_sequence_info2[,3]<-output_results_blast[,3]

######################################### This code block will be used to get the isolation source information for all sequeneces. In the file "Total_isolation.txt" we have the isolation source info from the database, and we are going to cluster that information based in out criteria

#isolation=as.data.frame(fread("Total_isolation.txt",header=F,sep="\t",na.strings=c("","NA")))
isolation=as.data.frame(fread("/storage/parras/databaseR/Tablas_taxa/Total_isolation.txt",header=F,sep="\t",na.strings=c("","NA")))

print(dim(isolation))
#isolation_criteria<-read.table("isolation_criteria.txt", na.strings=c("","NA"),row.names = 1,sep = '\t',header = F) #This file contains the clustering criteria into the different groups
isolation_criteria<-read.table("/home/parras/RESULTADOS_FINALES_NO_BORRAR/isolation_criteria.txt", na.strings=c("","NA"),row.names = 1,sep = '\t',header = F) #This file contains the clustering criteria into the different groups 

#isolation_exclusion_criteria<-read.table("isolation_exclusion_criteria.txt", na.strings=c("","NA"),row.names = 1,sep = '\t',header = F) #This file include an exclusion list for some sources, as they could be ambiguous. 
isolation_exclusion_criteria<-read.table("/home/parras/RESULTADOS_FINALES_NO_BORRAR/isolation_exclusion_criteria.txt", na.strings=c("","NA"),row.names = 1,sep = '\t',header = F) #This file include an exclusion list for some sources, as they could be ambiguous.  

isolation_criteria<-apply(isolation_criteria,2,tolower)
isolation_exclusion_criteria<-apply(isolation_exclusion_criteria,2,tolower)

isolation=isolation[isolation[,2] %in% sapply(strsplit(as.character(rownames(name_sequence_info2)), "\\."), "[[", 1),]

print(dim(isolation))

for(i in 1:length(name_sequence_info2[,1]))
{
  if(!isEmpty(isolation[isolation$V2==str_split(rownames(name_sequence_info2)[i],'\\.')[[1]][1],5]))
  {
    name_sequence_info2[i,4]<-isolation[isolation$V2==str_split(rownames(name_sequence_info2)[i],'\\.')[[1]][1],5][1]
  }
}

for(i in 1:length(output_results[,1]))
{
  output_results[i,13]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]],4],collapse = '?')
  output_results[i,14]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]],4],collapse = '?')
  
}

####We load isolation criteria data here

EnvironmentList <- list("Soil"=NULL,"Sediment"=NULL,"Plants&Fung"=NULL,"WasteWater"=NULL,"Build environment"=NULL,"Air"=NULL,"MarineW"=NULL,"Water"=NULL,"Aquaculture"=NULL,
"FreshW"=NULL,"Faecal"=NULL,"DomesticAnimals"=NULL,"Agricultur"=NULL,"Animals"=NULL,"Feed"=NULL,"Microbial"=NULL,"Misc"=NULL,"Human"=NULL,"Clinical"=NULL)
EnvironmentListS <- list("Soil_S"=NULL,"Sediment_S"=NULL,"Plants&Fung_S"=NULL,"WasteWater_S"=NULL,"Build environment_S"=NULL,"Air_S"=NULL,"MarineW_S"=NULL,"Water_S"=NULL,"Aquaculture_S"=NULL,
"FreshW_S"=NULL,"Faecal_S"=NULL,"DomesticAnimals_S"=NULL,"Agricultur_S"=NULL,"Animals_S"=NULL,"Feed_S"=NULL,"Microbial_S"=NULL,"Misc_S"=NULL,"Human_S"=NULL,"Clinical_S"=NULL)

EnvironmentList_Ex <- list("Soil_ex"=NULL,"Sediment_ex"=NULL,"Plants&Fung_ex"=NULL,"WasteWater_ex"=NULL,"Build environment_ex"=NULL,"Air_ex"=NULL,"MarineW_ex"=NULL,"Water_ex"=NULL,"Aquaculture_ex"=NULL,
"FreshW_ex"=NULL,"Faecal_ex"=NULL,"DomesticAnimals_ex"=NULL,"Agricultur_ex"=NULL,"Animals_ex"=NULL,"Feed_ex"=NULL,"Microbial_ex"=NULL,"Misc_ex"=NULL,"Human_ex"=NULL,"Clinical_ex"=NULL)
EnvironmentList_ExS <- list("Soil_exS"=NULL,"Sediment_exS"=NULL,"Plants&Fung_exS"=NULL,"WasteWater_exS"=NULL,"Build environment_exS"=NULL,"Air_exS"=NULL,"MarineW_exS"=NULL,"Water_exS"=NULL,"Aquaculture_exS"=NULL,
"FreshW_exS"=NULL,"Faecal_exS"=NULL,"DomesticAnimals_exS"=NULL,"Agricultur_exS"=NULL,"Animals_exS"=NULL,"Feed_exS"=NULL,"Microbial_exS"=NULL,"Misc_exS"=NULL,"Human_exS"=NULL,"Clinical_exS"=NULL)

for(i in 1:length(EnvironmentList))
{
	EnvironmentList[[i]]<-as.character(na.omit(as.character(as.matrix((isolation_criteria[rownames(isolation_criteria)==names(EnvironmentList)[i],])))))
	EnvironmentList_Ex[[i]]<-as.character(na.omit(as.character(as.matrix((isolation_exclusion_criteria[rownames(isolation_exclusion_criteria)==names(EnvironmentList)[i],])))))

	EnvironmentListS[[i]]<-EnvironmentList[[i]][nchar(EnvironmentList[[i]])< 5]
	EnvironmentList[[i]]<-EnvironmentList[[i]][nchar(EnvironmentList[[i]])> 4]
	
	EnvironmentList_ExS[[i]]<-EnvironmentList_Ex[[i]][nchar(EnvironmentList_Ex[[i]])< 5]
	EnvironmentList_Ex[[i]]<-EnvironmentList_Ex[[i]][nchar(EnvironmentList_Ex[[i]])> 4]

}

inclusion_clinical=c("calf","sponge","sole","fungal","fungus")


####Here we compute the number of matchers for each group

for(i in 1:length(name_sequence_info2[,1]))
{
  extracted_line<-tolower(name_sequence_info2[i,4])
 
	Matches_list<- list("matches_soil"=NULL,"matches_sediment"=NULL,"matches_plantsfungus"=NULL,"matches_wastewater"=NULL,"matches_buildenv"=NULL,"matches_air"=NULL,"matches_marinew"=NULL,
"matches_water"=NULL,"matches_aquaculture"=NULL,"matches_freshw"=NULL,"matches_faecal"=NULL,"matches_domesticanimals"=NULL,"matches_agricultur"=NULL,"matches_animals"=NULL,
"matches_feed"=NULL,"matches_microbial"=NULL,"matches_misc"=NULL,"matches_human"=NULL,"matches_clinical"=NULL)

	Matches_list_Ex<- list("matches_soil_ex"=NULL,"matches_sediment_ex"=NULL,"matches_plantsfungus_ex"=NULL,"matches_wastewater_ex"=NULL,"matches_buildenv_ex"=NULL,"matches_air_ex"=NULL,"matches_marinew_ex"=NULL,
"matches_water_ex"=NULL,"matches_aquaculture_ex"=NULL,"matches_freshw_ex"=NULL,"matches_faecal_ex"=NULL,"matches_domesticanimals_ex"=NULL,"matches_agricultur_ex"=NULL,"matches_animals_ex"=NULL,
"matches_feed_ex"=NULL,"matches_microbial_ex"=NULL,"matches_misc_ex"=NULL,"matches_human_ex"=NULL,"matches_clinical_ex"=NULL)

	for(o in 1:length(Matches_list))
	{
	  Matches_list[[o]]=grep(paste(EnvironmentList[[o]],collapse="|"),extracted_line, value=TRUE)
	  if(sort(str_split(extracted_line," ")[[1]] %in% EnvironmentListS[[o]],decreasing = TRUE)[1]&!is.na(extracted_line)){Matches_list[[o]]=extracted_line}
	  
	  Matches_list_Ex[[o]]=grep(paste(EnvironmentList_Ex[[o]],collapse="|"),extracted_line, value=TRUE)
	  if(sort(str_split(extracted_line," ")[[1]] %in% EnvironmentList_ExS[[o]],decreasing = TRUE)[1]&!is.na(extracted_line)){Matches_list_Ex[[o]]=extracted_line}
	  
	  if(length(Matches_list[[o]])>0 & length(Matches_list_Ex[[o]])<1)
	  {    
		Matches_list[[o]]<-names(EnvironmentList)[o]
	   }
	}  



if(length(Matches_list[[17]])>0 &
     length(Matches_list[[1]])<1 & 
      length(Matches_list[[2]])<1 &   
      length(Matches_list[[3]])<1  &
     length(Matches_list[[4]])<1  &
     length(Matches_list[[5]])<1   &
     length(Matches_list[[6]])<1 &  
     length(Matches_list[[7]])<1   &
     length(Matches_list[[8]])<1 &  
     length(Matches_list[[9]])<1 &  
     length(Matches_list[[10]])<1  & 
     length(Matches_list[[11]])<1   &
     length(Matches_list[[12]])<1 & 
     length(Matches_list[[13]])<1 &  
     length(Matches_list[[14]])<1 &
     length(Matches_list[[15]])<1  &
     length(Matches_list[[16]])<1 & 
     length(Matches_list[[17]])<1& 
     length(Matches_list[[18]])<1 & 
     length(Matches_list_Ex[[17]])<1)
  {
	Matches_list[[17]]<-"Misc"
  }
  
  if(length(Matches_list[[18]])>0  & length(Matches_list[[14]])<1 & length(Matches_list_Ex[[18]])<1)
  {
	Matches_list[[18]]<-"Human"
    
  }
  
  if(length(Matches_list[[19]])>0 & length(Matches_list[[14]])<1 & length(Matches_list[[3]])<1 & length(Matches_list_Ex[[19]])<1)
  {
	Matches_list[[19]]<-"Clinical"  
  }else if(length(Matches_list[[19]])>0 & (sort(str_split(extracted_line," ")[[1]] %in% inclusion_clinical,decreasing = TRUE)[1]&!is.na(extracted_line)))
  {
    Matches_list[[19]]<-"Clinical"
  }
  
  Matches_list[Matches_list=="character(0)"]=NULL
  name_sequence_info2[i,5]<-paste3(Matches_list,collapse = '?') #Here we save the results for each sequences, which isolation source clusters are any of these found						 
}

for(i in 1:length(output_results[,1])) #We get the results into the axuliary table where we have all the results together
{
  extracted_line<-str_split(output_results[i,13],"\\?")[[1]]
  isolation_output_table=NULL
  
  for(o in 1:length(extracted_line))
  {
    isolation_output_table=c(isolation_output_table,name_sequence_info2[name_sequence_info2$V4 %in% extracted_line[o],5][1])
  }
  
  
  #holi=data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE)))
  #holi=subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))
  output_results[i,15]=paste(subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))$Var1,subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))$Freq,collapse = ";")

  extracted_line<-str_split(output_results[i,14],"\\?")[[1]]
  isolation_output_table=NULL
  
  for(o in 1:length(extracted_line))
  {
    isolation_output_table=c(isolation_output_table,name_sequence_info2[name_sequence_info2$V4 %in% extracted_line[o],5][1])
  }
  
  #holi=data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE)))
  #holi=subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))
  output_results[i,16]=paste(subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))$Var1,subset(data.frame(table(gsub(" ","",unlist(str_split(paste(isolation_output_table,sep = ""),pattern = ",")),ignore.case = TRUE))), grepl('[aA-zZ]', Var1) & !grepl('NA', Var1))$Freq,collapse = ";")
}
#################################This block is to compute pathogenic bacteria info

name_sequence_info2[,6]<-NA
name_sequence_info2[,7]<-NA
output_results[i,17]<-NA
output_results[i,18]<-NA

################### Get the info from Blat results into the results table

output_results[,19:22]<-NA

for(i in 1:length(output_results[,1]))
{
  if(!isEmpty(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]>90,2])))
  {
    output_results[i,19]<-paste3(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]>90,2]),collapse = '?')
  }
  if(!isEmpty(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]<90,2])))
  {
    output_results[i,21]<-paste3(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]<90,2]),collapse = '?')
  }
  if(!isEmpty(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]>90,3]))
  {
    output_results[i,20]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]>90,3],collapse = '?')
  }
  if(!isEmpty(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]<90,3]))
  { 
    output_results[i,22]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,4]),";")[[1]] & name_sequence_info2[,3]<90,3],collapse = '?')
  }
}

output_results[,23:26]<-NA

for(i in 1:length(output_results[,1]))
{
  if(!isEmpty(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]>90,2])))
  {
    output_results[i,23]<-paste3(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]>90,2]),collapse = '?')
  }
  if(!isEmpty(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]<90,2])))
  {
    output_results[i,25]<-paste3(as.character(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]<90,2]),collapse = '?')
  }
  if(!isEmpty(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]>90,3]))
  {
    output_results[i,24]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]>90,3],collapse = '?')
  }
  if(!isEmpty(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]<90,3]))
  { 
    output_results[i,26]<-paste3(name_sequence_info2[rownames(name_sequence_info2) %in% str_split(as.character(output_results[i,5]),";")[[1]] & name_sequence_info2[,3]<90,3],collapse = '?')
  }
}

########## Save the distance results into the table

pairwise_table<-output_results_pairwise

output_results[,27]<-pairwise_table[,2]
output_results[,28]<-pairwise_table[,3]

########## Now we process all the information to get a tree plotted with all the information we have about each sequence

phyla<-c("Abditibacteriota","Acidobacteria","Actinobacteria","Armatimonadetes","Bacteroidetes",
         "Balneolaeota","Caldiserica","Calditrichaeota","Candidate division KSB1","Candidate division TM6",
         "Candidate division WOR-1","Candidate division WPS-2","Candidatus Accumulibacter","Candidatus Aminicenantes","Candidatus Competibacteraceae",
         "Candidatus Cryosericota","Candidatus Delongbacteria","Candidatus Dependentiae","Candidatus Eisenbacteria","Candidatus Eremiobacteraeota",
         "Candidatus Kapabacteria","Candidatus Kerfeldbacteria","Candidatus Kryptonia","Candidatus Latescibacteria","Candidatus Lindowbacteria",
         "Candidatus Margulisbacteria","Candidatus Marinimicrobia","Candidatus Melainabacteria","Candidatus Omnitrophica","Candidatus Peregrinibacteria",
         "Candidatus Saccharibacteria","Candidatus Sumerlaeota","Candidatus Wallbacteria","Chlamydiae","Chlorobi",
         "Chloroflexi","Cyanobacteria","Deinococcus-Thermus","Elusimicrobia","Fibrobacteres",
         "Firmicutes","Fusobacteria","Gemmatimonadetes","Ignavibacteriae","Lentisphaerae",
         "Nitrospirae","Nitrospinae","Planctomycetes","Proteobacteria","Rhodothermaeota","Spirochaetes",
         "Synergistetes","Tenericutes","Verrucomicrobia","Thermotogae","Dictyoglomi")

colours_phyla<-c("#FCFFA4","#991914","#267C03","#4c4cb2","#E2AC6F",
                 "#FFA6AF","#9eff99","#8d8af2","#000000","#000000",
                 "#000000","#000000","#000000","#000000","#000000",
                 "#000000","#000000","#000000","#000000","#000000",
                 "#000000","#000000","#000000","#000000","#000000",
                 "#000000","#000000","#000000","#000000","#000000",
                 "#000000","#000000","#000000","#E0FF99","#89f3ff",
                 "#f49b0c","#f4b2bd","#30840e","#e81078","#f2ca8a",
                 "#6C74AF","#33ce5a","#db677a","#f458da","#2667ff",
                 "#ff7563","#B983FF","#7bed82","#E83B68","#eec9ff","#A5C1E5",
                 "#d3442e","#574cce","#DDD547","#8C543A","#8B3E2F")

color_summed<-as.list(setNames(colours_phyla,phyla))

name_fil<-data.frame(name_sequence_info2[,1])
row.names(name_fil)<-row.names(name_sequence_info2)
colnames(name_fil)<-"Phylum"

for(w in 1:length(unique(name_fil[,1])))
{
	if(isEmpty(grep(".",color_summed[as.character(unique(name_fil[,1])[w])][[1]])) & !is.na(as.character(unique(name_fil[,1])[w])))
	{
		color_summed[as.character(unique(name_fil[,1])[w])]="#000000"
	}
}

name_per<-data.frame(name_sequence_info2[,3])
row.names(name_per)<-row.names(name_sequence_info2)
colnames(name_per)<-"Blast best hit %"
name_prot<-data.frame(name_sequence_info2[,2])
row.names(name_prot)<-row.names(name_sequence_info2)
name_prot[,1]=as.character(name_prot[,1])

for(i in 1:length(name_prot[,1]))
{
  if(!isEmpty(grep("-",as.character(name_prot[i,1]))))
  {
    name_prot[i,1]=sapply(strsplit(as.character(name_sequence_info2[i,2]), "-"), "[[", 1)
  }else
  {
    name_prot[i,1]=sapply(strsplit(as.character(name_sequence_info2[i,2]), "_"), "[[", 1)
  }
}								 
colnames(name_prot)<-"Protein"
name_source<-data.frame(name_sequence_info2[,5])
row.names(name_source)<-row.names(name_sequence_info2)
colnames(name_source)<-"Isolation source"

tree_colors<-data.frame(matrix(ncol = 5,nrow = 26))
rownames(tree_colors)<-c("class_A.hmm","class_B_1_2.hmm","class_C.hmm","class_B_3.hmm","class_D_1.hmm","class_D_2.hmm","qnr.hmm","tet_efflux.hmm","tet_enzyme.hmm",
"tet_rpg.hmm","mph.hmm","erm_typeF.hmm","B1.hmm","class_d_1_2","methyltransferase_grp_1_2",
"aac2p.hmm","aac3_class1.hmm","aac3_class2.hmm","aac6p_class1.hmm","aac6p_class2.hmm","aac6p_class3.hmm","aph2b.hmm","aph3p.hmm","aph6.hmm","aac6p_complete","16S_RMT")

tree_colors[,5]=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
tree_colors[,1]=c(1,1,1,1,1,1,1,1,0.5,1,1,1,1,1,1,1,1,0.5,1,1,1,1,1,1,1,0.5)
tree_colors[,2]=c(16,10,14.5,18.5,11,11,11,16,8.5,30,14,11,11,33.5,12.5,5.5,5.5,11,15,19,6,10.5,14,20,25,5)
tree_colors[,3]=c(23.5,14,20,26,16,16,15.5,23,12,42.5,19.5,16,16,47,17,7.5,7.5,15.5,21,27,8,14.5,20,28,35.5,7.5)
tree_colors[,4]=c(140,65,105,160,65,235,80,145,18,215,100,65,85,205,160,150,150,150,150,150,150,150,150,150,150,150)

elemento_impr=model_list[j]

p8<-ggtree(tree_l,layout = 'circular',branch.length='none') +
  geom_tiplab(size=0) + geom_point2(aes(subset=(node %in% Nodes_to_color)),color="green",size=2) + geom_nodelab2(aes(subset=(node %in% Nodes_to_color)),color="red",size=5)

p8 <- gheatmap(p8, name_fil, offset=tree_colors[row.names(tree_colors)==elemento_impr,1], width=.1,
               colnames_angle=90, colnames_offset_y = .25,color = NULL) +
 scale_fill_manual(values=color_summed)


p8 <- p8 + new_scale_fill()
p8<-gheatmap(p8, name_per, offset=tree_colors[row.names(tree_colors)==elemento_impr,3], width=.05,
             colnames_angle=90, colnames_offset_y = .25,color=NULL) +
  scale_fill_viridis_c(option="B", name="Blast best hit %")

p8<- p8 + new_scale_fill()
p8<-gheatmap(p8, name_prot, offset=tree_colors[row.names(tree_colors)==elemento_impr,2], width=.05,
             colnames_angle=90, colnames_offset_y = .25,color = NULL) +
  scale_fill_manual(values=randomColor(length(levels(factor(name_prot[,1])))))  #+ theme(legend.position = "none") #####DEJAR DE SILENCIAR SI HAY MUCHAS PROT




pdf("Circular_tree_plot.pdf",height = 20, width = 20) #All information is printed and saved
plot(p8)
dev.off()

name_sequence_info3<-name_sequence_info2
colnames(name_sequence_info3)=c("Phylum","Related protein","% Similarity","Isolation","Isolation_cluster","Pathogenicity","MGE")
write.table(name_sequence_info3,"Metadata_results.txt",sep = "\t",quote = F)

output_results2<-output_results
colnames(output_results2)=c("Node","Phylum1","Phylum2","Species1","Species2","NSpe1","NSpe2","Block1","Block2","Block3","NSeqsPlasmids1","NSeqsPlasmids2",
                        "IsolationSource1","IsolationSource2","Conteo1","Conteo2","Patho1","Patho2","Known Prot 1","Identity - known prot 1","New Prot 1",
                        "Identity - new prot 1","Known Prot 2","Identity - known prot 2","New Prot 2","Identity - new prot 2","Min Pair","Max Pair")
write.table(output_results2,"Events_output.txt",sep = "\t",row.names = F,quote = F)

output_results_min<-output_results2[,c(1,2,3,6,7,11,12,13,14,15,16,17,18,20,22,24,26,27,28)]
write.table(output_results_min,"Events_output_min.txt",sep = "\t",row.names = F,quote = F)

#system(paste("mv Blast_Antibiotic_analysis Taxonomy_info_dataset.txt jolin.txt prueba_a_ver_ids_mod.txt Alignment_file resultados_discrepacias_isolation_min3.txt resultados_discrepacias_isolation3.txt resultados_prot3.txt Resistance_genes_list_species.fasta Resistance_genes_list.fa Tree_file tree_paint_jj_rect_filo_pato3.pdf ./",model_list[j],"/",sep=""))

system(paste("mv Blast_Antibiotic_analysis Taxonomy_info_dataset.txt Alignment_file Events_output_min.txt Events_output.txt Metadata_results.txt Resistance_genes_list_species.fasta Resistance_genes_list.fa Tree_file Circular_tree_plot.pdf ./",model_list[j],"/",sep=""))
}