library(shiny)
library(tidyverse)
library(ggbeeswarm)

shiny_log <- "shiny_responses.tsv"
fields <- c("mutation","Rating","User","tag","comment")
review_results = suppressMessages(read_tsv(shiny_log,col_names=fields)) 


fix_alt = function(Ref,Alt,kept_alt,basename){
  case_when(
    kept_alt == "snapshot.png" & str_detect(basename,paste0(Ref,"----")) ~ "-", #deletion
    kept_alt == "snapshot.png" & str_detect(basename,paste0("----",Ref)) ~ Ref, #insertion
    TRUE ~ Alt
  )
}

fix_ref = function(Ref,Alt,kept_alt,basename){
  case_when(
    kept_alt == "snapshot.png" & str_detect(basename,paste0(Ref,"----")) ~ Ref, #deletion
    kept_alt == "snapshot.png" & str_detect(basename,paste0("----",Ref)) ~ "-", #insertion
    TRUE ~ Ref
  )
}

# load the file names for screenshots of variants unique to Reddy or GAMBL (mostly)
options_df = data.frame(full=dir(recursive=T,pattern=".png")) %>% 
  dplyr::filter(!grepl("pairs",full)) %>%
  mutate(base=basename(full)) %>%
  mutate(basename=base) %>% 
  mutate(base=str_replace(base,"HLA-","HLA_")) %>%
  mutate(base=str_replace(base,"MEF2BNB-","MEF2BNB_")) %>%
  separate(base,into=c("Region","End","Gene","sample_id","Ref","Alt"),sep="-+") %>% 
  mutate(chrpos=Region) %>%
  separate(chrpos,into=c("Chromosome","Start"),sep=":") %>%
  mutate(Start=as.numeric(Start)) %>%
  mutate(End = as.numeric(End)) %>%
  mutate(gap=End-Start+1) %>%
  mutate(Start_Position = case_when(gap ==402 ~ Start + 200,
                                    gap > 150 ~ Start + 75,
                                    TRUE ~ Start)) %>% 
  mutate(kept_alt=Alt) %>%
  #mutate(Ref=ifelse(grepl(paste0(kept_alt,"----"),basename),Ref,"-")) %>%
  mutate(Alt=fix_alt(Ref,Alt,kept_alt,basename)) %>%
  mutate(Ref=fix_ref(Ref,Alt,kept_alt,basename)) %>%
  select(-kept_alt) %>%
  arrange(Gene,Region) 


# load the file names for screenshots of variants shared by Reddy and GAMBL (mostly)
#process slightly differnetly due to a change in the file naming 
paired_df = data.frame(full=dir(recursive=T,pattern=".png")) %>% 
  dplyr::filter(grepl("pairs",full)) %>%
  mutate(base=basename(full)) %>%
  mutate(basename=base) %>% 
  mutate(base=str_replace(base,"HLA-","HLA_")) %>%
  mutate(base=str_replace(base,"MEF2BNB-","MEF2BNB_")) %>%
  mutate(base=str_remove(base,".pairs")) %>%
  separate(base,into=c("Region","End","Gene","sample_id","Ref","Alt"),sep="-+") %>% 
  mutate(chrpos=Region) %>%
  separate(chrpos,into=c("Chromosome","Start"),sep=":") %>%
  mutate(Start=as.numeric(Start)+1) %>%
  mutate(End = as.numeric(End)) %>%
  #mutate(Start_Position = Start) %>%
  mutate(gap=End-Start+1) %>%
  mutate(Start_Position = case_when(gap ==402 ~ Start + 200,
                                    gap > 150 ~ Start + 75,
                                    TRUE ~ Start)) %>% 
  arrange(Gene,Region)

options_df = bind_rows(paired_df,options_df) %>% arrange(Gene,Region)


# load all the variant reviews so they can be matched up to the original variant calls
review_all = full_join(options_df,review_results,by=c("basename"="mutation")) %>% 
  select(-full,-Region,-End,-Start,-gap,-comment,-tag) %>% 
  unique() %>%
  group_by(basename) %>%
  mutate(Mean_Rating=mean(Rating)) %>%
  mutate(All_Rating=paste(Rating,collapse=",")) %>%
  slice_head() %>% 
  ungroup() %>% 
  select(-User,-basename) %>%
  filter(!Gene %in% c("MEF2BNB_MEF2B")) %>% 
  mutate(Chromosome = str_remove(Chromosome,pattern = "chr"))
    


#match up calls reported in Reddy study so they can be handled on their own
#This file is just the SNV
reddy_annotated = read_tsv("data/reddy_snv_all_reannotated.maf.gz")


#some are missing from this file. Use the one above instead.
reddy_orig = read_tsv("data/mutations_reddy_original_annotations_reformatted.maf.gz") %>% 
  mutate(Tumor_Sample_Barcode = paste0("Reddy_",Tumor_Sample_Barcode,"T")) %>% 
  mutate(Chromosome = str_remove(Chromosome,"chr"))

reddy_indel = filter(reddy_orig,Variant_Type != "SNP") %>% 
  filter(!Hugo_Symbol %in% c("MEF2BNB-MEF2B")) %>% unique()


annotated_snv_from_reddy = left_join(select(review_all,-Gene),
                                 reddy_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
                                                      "Chromosome"="Chromosome",
                                                      "Start_Position"="Start_Position")) %>%
  filter(!is.na(Hugo_Symbol))


annotated_indel_from_reddy = left_join(select(review_all,-Gene),
                                       reddy_indel,by=c("sample_id"="Tumor_Sample_Barcode",
                                                          "Chromosome"="Chromosome",
                                                          "Start_Position"="Start_Position")) %>%
  filter(!is.na(Hugo_Symbol)) %>% 
  
  unique()



#write_tsv(select(review_all,Chromosome,Start_Position,Rating,User),file="all_reviews_Reddy.grch37.tsv")

#match up calls from GAMBL so they can be handled on their own

gambl_annotated = read_tsv("data/reddy_mutations_gambl_dlbclgenes.maf.gz") 
gambl_annotated = group_by(gambl_annotated,Chromosome,Start_Position,Tumor_Sample_Barcode) %>% 
  slice_head(n=1) %>% ungroup()

#exclude indels since they aren't being matched properly anyway
gambl_annotated = dplyr::filter(gambl_annotated,Variant_Type %in% c("SNP","DNP","TNP"))


#this is the intersect
reddy_annotated_in_gambl = left_join(dplyr::select(gambl_annotated,
                                                   Chromosome,Start_Position,Tumor_Sample_Barcode, Tumor_Seq_Allele2),
                                     reddy_annotated) %>%
  dplyr::filter(!is.na(Hugo_Symbol))


annotated_intersect = inner_join(review_all,
                                dplyr::select(reddy_annotated_in_gambl,-Gene),by=c("sample_id"="Tumor_Sample_Barcode",
                                                              "Chromosome"="Chromosome",
                                                              "Start_Position"="Start_Position")) 

#annotated_from_gambl = left_join(select(review_all,-Gene),
#                                 gambl_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
#                                                      "Chromosome"="Chromosome",
#                                                      "Start_Position"="Start_Position")) %>%
#  filter(!is.na(Hugo_Symbol),!is.na(User))


#GAMBL variants not also covered by Reddy
unannotated_from_gambl = anti_join(gambl_annotated,annotated_intersect,by=c("Tumor_Sample_Barcode"="sample_id",
                                                                          "Chromosome"="Chromosome",
                                                                          "Start_Position"="Start_Position"))
annotated_from_gambl = inner_join(review_all,unannotated_from_gambl,by=c("sample_id"="Tumor_Sample_Barcode",
                                                                         "Chromosome"="Chromosome",
                                                                         "Start_Position"="Start_Position"))
all_annotated_from_gambl = inner_join(review_all,gambl_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
                                                                         "Chromosome"="Chromosome",
                                                                         "Start_Position"="Start_Position"))
reddy_g = dplyr::filter(lymphoma_genes_dlbcl_v_latest,
                        Reddy==TRUE) %>% pull(Gene)
all_annotated_from_gambl_reddygenes = dplyr::filter(all_annotated_from_gambl,
                                                    Hugo_Symbol %in% reddy_g)
#This still includes the intersect. Need to remove the contents of annotated_intersect

gambl_only_reddy_genes = anti_join(all_annotated_from_gambl_reddygenes,annotated_intersect,by=c("sample_id",
                                                                                                "Start_Position",
                                                                                                "Alt"))

reddy_only_snv = anti_join(annotated_snv_from_reddy,annotated_intersect,by=c("sample_id",
                                                                "Chromosome",
                                                                "Start_Position"))

reddy_only_indel = anti_join(annotated_indel_from_reddy,annotated_intersect,by=c("sample_id",
                                                                             "Chromosome",
                                                                             "Start_Position"))


unreviewed = anti_join(gambl_annotated,review_all,by=c("Tumor_Sample_Barcode"="sample_id",
                                                       "Chromosome"="Chromosome",
                                                       "Start_Position"="Start_Position"))

write_tsv(select(unreviewed,Hugo_Symbol,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode),file="unreviewed_gambl.maf")

#SNV comparison
reddy_only_snv = mutate(reddy_only_snv,group="Reddy_only")
all_annotated_from_gambl = mutate(all_annotated_from_gambl,group="GAMBL_only")
reddy_v_gambl_snv_rating = bind_rows(select(reddy_only_snv,Chromosome,Start_Position,sample_id,Mean_Rating,group),
                                     select(all_annotated_from_gambl,Chromosome,Start_Position,sample_id,Mean_Rating,group))
annotated_intersect = mutate(annotated_intersect,group="Intersect")
threeway_compare = bind_rows(reddy_v_gambl_snv_rating,
                             select(annotated_intersect,Start_Position,sample_id,Mean_Rating,group))
#ggplot(reddy_v_gambl_snv_rating,aes(x=Mean_Rating)) + geom_histogram() + facet_wrap(~group,ncol=1,scales="free_y") + theme_cowplot()
#ggsave("rating_SNV_comparison_Reddy_vs_GAMBL_no_intersect.pdf")
#save them for downstream analysis

ggplot(threeway_compare,aes(x=Mean_Rating)) + geom_histogram() + facet_wrap(~group,ncol=1,scales="free_y") + theme_cowplot()
ggsave("rating_SNV_comparison_Reddy_vs_GAMBL_with_intersect.pdf")


dplyr::select(all_annotated_from_gambl_reddygenes,-Gene,-Ref,-Alt,-Rating,-coding_mutations,-MeanCorrectedCoverage) %>%
  dplyr::rename("Tumor_Sample_Barcode"="sample_id") %>%
  write_tsv(file="gambl_only_mutations_in_Reddy_genes_with_review.maf")

dplyr::select(all_annotated_from_gambl,-Gene,-Ref,-Alt,-Rating,-coding_mutations,-MeanCorrectedCoverage) %>%
  dplyr::rename("Tumor_Sample_Barcode"="sample_id") %>%
  write_tsv(file="gambl_only_mutations_lymphoma_genes_with_review.maf")


dplyr::select(reddy_only_snv,-Ref,-Alt) %>%
  dplyr::rename("Tumor_Sample_Barcode"="sample_id") %>%
  write_tsv(file="reddy_only_SNVs_with_review.maf")


#dplyr::select(annotated_intersect) %>%
  #dplyr::rename("Tumor_Sample_Barcode"="sample_id") %>%
write_tsv(annotated_intersect,file="reddy_gambl_shared_mutations_with_review.maf")


minimal_rating_info = bind_rows(dplyr::select(annotated_from_gambl,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,sample_id,Rating,tag,comment),
          dplyr::select(annotated_intersect,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,sample_id,Rating,tag,comment),
          dplyr::select(annotated_from_reddy,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,sample_id,Rating,tag,comment)) %>% unique()

write_tsv(minimal_rating_info,file="manualreview--grch37--capture.tsv")

annotated_from_gambl_score = dplyr::select(annotated_from_gambl,Hugo_Symbol,Rating) %>% mutate(group="GAMBL-only")
annotated_from_reddy_score = dplyr::select(annotated_from_reddy,Hugo_Symbol,Rating) %>% mutate(group="Reddy-only")
intersect_score = dplyr::select(annotated_intersect,Hugo_Symbol,Rating) %>% mutate(group="GAMBL_and_Reddy")

compare_three = bind_rows(annotated_from_gambl_score,annotated_from_reddy_score,intersect_score) %>% 
  arrange(Hugo_Symbol)

ui <- fluidPage(
  titlePanel("Mutation review status"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Choose a Gene", choices = c("All",unique(compare_three$Hugo_Symbol)))
    ),
    mainPanel(
      #tableOutput("all_reviews"),
      tags$style(
        ".calc1{
        color:red;
        position:absolute;
        margin-left: -220px;
        margin-top: -200px
        }
        .calc2{
        color:blue;
        position: absolute;
        margin-left: -220px;
        margin-top: -215px
        }
        "
      ),
      div(plotOutput("graph1"))
      #div(class="calc1", textOutput("calc1")),
      #br(),
      #div(plotOutput("graph2")),
      #div(class="calc2", textOutput("calc2"))
    )
  )
)
server <- function(input, output, session) {
  g = reactive({
    if(input$gene=="All"){
      compare_three
    }else{
      compare_three %>% dplyr::filter(Hugo_Symbol==input$gene)
    }
  })
  
  
  output$graph1<-renderPlot({
    #zero value placeholder for filling in missing bins
    blank_df = data.frame(Rating=c(0,1,2,3,4,5),n=c(0,0,0,0,0,0))
    
    #actual counts for bins with at least one value (some may be missing)
    g() %>% ggplot(aes(y=Rating)) + 
      geom_histogram() + 
      facet_wrap(~group) 
  })
  #output$calc1 <- renderText({ 
    #actual counts for bins with at least one value (some may be missing)
  #  counted_df = r() %>% 
  #    group_by(Rating) %>%
  #    tally()
    
    #Calculating mean and median
  #  med_data = rep(counted_df$Rating, counted_df$n)
  #  mean_df = mean(med_data)
  #  median_df = median(med_data)
    
  #  paste("Mean : ", round(mean_df, digits = 3),"Median : ", median_df)
    
  #})

  
}
shinyApp(ui, server)
