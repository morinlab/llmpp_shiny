library(shiny)
library(tidyverse)
library(ggbeeswarm)

shiny_log <- "shiny_responses.tsv"
fields <- c("mutation","Rating","User","tag","comment")
review_results = suppressMessages(read_tsv(shiny_log,col_names=fields)) 


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
  arrange(Gene,Region) 

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
  arrange(Gene,Region)

options_df = bind_rows(paired_df,options_df) %>% arrange(Gene,Region)



review_all = full_join(options_df,review_results,by=c("basename"="mutation")) %>% 
  mutate(Start=as.numeric(Start),End=as.numeric(End)) %>%
  mutate(Start_Position = Start + 200) %>% 
  group_by(basename) %>%
  mutate(Rating=mean(Rating)) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  mutate(Chromosome = str_remove(Chromosome,pattern = "chr"))
#View(review_results)
#match up calls reported in Reddy study so they can be handled on their own
reddy_annotated = read_tsv("data/reddy_snv_all_reannotated.maf.gz")




annotated_from_reddy = left_join(select(review_all,-Gene),
                                 reddy_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
                                                      "Chromosome"="Chromosome",
                                                      "Start_Position"="Start_Position")) %>%
  filter(!is.na(Hugo_Symbol),!is.na(User))



#match up calls from GAMBL so they can be handled on their own

gambl_annotated = read_tsv("data/reddy_mutations_gambl_dlbclgenes.maf.gz") 

#exclude indels since they aren't being matched properly anyway
gambl_annotated = dplyr::filter(gambl_annotated,Variant_Type %in% c("SNP","DNP","TNP"))



reddy_annotated_in_gambl = left_join(dplyr::select(gambl_annotated,Chromosome,Start_Position,Tumor_Sample_Barcode, Tumor_Seq_Allele2),reddy_annotated) %>%
  dplyr::filter(!is.na(Hugo_Symbol))


annotated_intersect = inner_join(review_all,
                                dplyr::select(reddy_annotated_in_gambl,-Gene),by=c("sample_id"="Tumor_Sample_Barcode",
                                                              "Chromosome"="Chromosome",
                                                              "Start"="Start_Position")) 
#%>%
#  filter(!is.na(Hugo_Symbol),!is.na(User))

annotated_from_gambl = left_join(select(review_all,-Gene),
                                 gambl_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
                                                      "Chromosome"="Chromosome",
                                                      "Start_Position"="Start_Position")) %>%
  filter(!is.na(Hugo_Symbol),!is.na(User))

dplyr::select(annotated_from_gambl,-User,-Start,-End,-Ref,-Alt,-full,-Region,-basename) %>%
  rename("sample_id"="Tumor_Sample_Barcode") %>%
  write_tsv(file="gambl_only_mutations_with_review.maf")


dplyr::select(annotated_from_reddy,-User,-Start,-End,-Ref,-Alt,-full,-Region,-basename) %>%
  rename("sample_id"="Tumor_Sample_Barcode") %>%
  write_tsv(file="reddy_only_mutations_with_review.maf")


dplyr::select(annotated_intersect,-User,-Start,-End,-Ref,-Alt,-full,-Region,-basename) %>%
  rename("sample_id"="Tumor_Sample_Barcode") %>%
  write_tsv(file="reddy_gambl_shared_mutations_with_review.maf")

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
