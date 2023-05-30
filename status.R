library(shiny)
library(tidyverse)

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


review_all = full_join(options_df,review_results,by=c("basename"="mutation")) %>% 
  mutate(Start=as.numeric(Start),End=as.numeric(End)) %>%
  mutate(Start_Position = Start + 200) 

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

annotated_from_gambl = left_join(select(review_all,-Gene),
                                 gambl_annotated,by=c("sample_id"="Tumor_Sample_Barcode",
                                                      "Chromosome"="Chromosome",
                                                      "Start_Position"="Start_Position")) %>%
  filter(!is.na(Hugo_Symbol),!is.na(User))

#gene_id = "MYD88"
#filter(annotated_from_reddy,Hugo_Symbol==gene_id) %>% 
#  group_by(Rating) %>% tally() %>%
#  ggplot(aes(x=Rating,y=n)) + geom_col() + xlim(c(0,5)) + ggtitle(gene_id)

ui <- fluidPage(
  titlePanel("Mutation review status"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Choose a Gene", choices = c("All",unique(annotated_from_reddy$Hugo_Symbol)))
    ),
    mainPanel(
      #tableOutput("all_reviews"),
      plotOutput("graph1"),
      plotOutput("graph2")
    )
  )
)
server <- function(input, output, session) {
  g = reactive({
    if(input$gene=="All"){
      annotated_from_gambl
    }else{
      annotated_from_gambl %>% dplyr::filter(Hugo_Symbol==input$gene)
    }
  })
  r = reactive({
    if(input$gene=="All"){
      annotated_from_reddy
    }else{
      annotated_from_reddy %>% dplyr::filter(Hugo_Symbol==input$gene)
    }
  })
  

  output$graph1<-renderPlot({
    #zero value placeholder for filling in missing bins
    blank_df = data.frame(Rating=c(0,1,2,3,4,5),n=c(0,0,0,0,0,0))
    
    #actual counts for bins with at least one value (some may be missing)
    counted_df = r() %>% 
      group_by(Rating) %>%
      tally() 
    
    #bind the two and take the highest value per bin
    bind_rows(blank_df,counted_df) %>% 
      group_by(Rating) %>% 
      arrange(desc(n)) %>% 
      slice_head() %>%
      ggplot(aes(x=Rating,y=n)) + 
      geom_col() + 
      ggtitle(paste(input$gene, "SNVs from Reddy"))
  })
  
  output$graph2<-renderPlot({
    #zero value placeholder for filling in missing bins
    blank_df = data.frame(Rating=c(0,1,2,3,4,5),n=c(0,0,0,0,0,0))
    
    #actual counts for bins with at least one value (some may be missing)
    counted_df = g() %>% 
      group_by(Rating) %>%
      tally() 
    
    #bind the two and take the highest value per bin
    bind_rows(blank_df,counted_df) %>% 
      group_by(Rating) %>% 
      arrange(desc(n)) %>% 
      slice_head() %>%
      ggplot(aes(x=Rating,y=n)) + 
      geom_col() + 
      ggtitle(paste(input$gene, "SNVs from GAMBL"))
  })
}
shinyApp(ui, server)

