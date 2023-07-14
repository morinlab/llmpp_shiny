library(shiny)
library(tidyverse)
library(shinyjs)

# Define the fields we want to save from the form
fields <- c("mutation","rb","user","tag","comment")

shiny_log <- "shiny_responses.tsv"

#try to get the git user ID
get_user = function(){
  # on the server automatically set a user ID that should be changed if users want their reviews logged
  return("anonymous")
}


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



save_data=function(data){
  if(class(data[["tag"]]) == "character"){
    data[["tag"]] = unlist(paste0(data[["tag"]],collapse=";"))
  }else{
    data[["tag"]] = ""
  }
  
  data <- data.frame(t(sapply(data,c)))
  colnames(data)=fields
  
  suppressMessages(write_tsv(data,file=shiny_log,append = T))
  review_results = suppressMessages(read_tsv(shiny_log,col_names=fields))
  return(review_results)
}

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Welcome, Mutation Reviewer!"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Choose a Gene or hit the button below", choices = unique(options_df$Gene)),
      actionButton("random", "Suggest a random gene"),
      radioButtons("status","Which group of variants do you want to review?", choiceNames=list("Unreviewed","Reviewed"),choiceValues=list("unreviewed","reviewed")),
      selectInput("mutation", "Pick a Mutation", choices = NULL),
      radioButtons("rb", "Mutation quality:",width="100%",
                   choiceNames = list(
                     "0: Zero support for the variant in the reads",
                     "1: Minimal support and/or severe confounders",
                     "2: Low support (2-3 molecules) or other confounders",
                     "3: Modest support but some uncertainty or a confounder",
                     "4: Good support",
                     "5: Excellent support, no ambiguity"
                   ),
                   choiceValues = list("0", "1", "2","3","4","5")
      ),
      checkboxGroupInput(
        "tag",
        "Select any tags that apply (Optional)",
        choiceNames = c("Adjacent indel: attributable to misalignment caused by a nearby insertion or deletion",
                        "Directional: Variant is only/mostly found on reads in the same direction",
                        "Multiple Variants: Variant appears multi-allelic, having read support for two or more non-reference alleles",
                        "Low Variant Frequency: Low VAF in tumour",
                        "End of reads: Variant is only seen within ~30 bp of end of all variant-supporting reads",
                        "High Discrepancy Region (supported by reads that have other recurrent mismatches across the track)",
                        "Multiple Mismatches: Variant is supported by reads that have other mismatched base pairs",
                        "Low Mapping quality: Variant is mostly supported by reads that have low mapping quality",
                        "Short Inserts Only: Variant is exclusively found on short fragments such that sequencing from each end results in overlapping reads",
                        "Low Count Tumor: Variant has inadequate coverage in the tumour",
                        "Same Start End: Variant is only observed in reads that start and stop at the same positions",
                        "No Count Normal: No coverage in the matched normal or data from matched normal is unavailable",
                        "Mononucleotide repeat: adjacent to a region in the genome that has a single-nucleotide repeat (e.g., AAAAAA…)",
                        "Dinucleotide repeat: adjacent to a region in the genome that has two alternating nucleotides (e.g. TGTG...)",
                        "Tandem repeat: adjacent to a region in the genome that has three or more alternating nucleotides (e.g., GTGGTGGTG…)",
                        "Ambiguous other: surrounded by inconclusive genomic features that cannot be explained by other tags"),
        choiceValues = c("AI","D","MV","LVF","E","HDR","MM","LM","SIO","LCT","SSE",
                         "NCN","MN","DN","TR",
                         "AO"),
        selected = c("NCN")
        ,inline = F
      ),
      tableOutput("leaderboard")
    ),
    mainPanel(
      textInput("comment","Enter a comment (optional)"),
      
      textInput("user", "Please change this to your GitHub user ID:",value = get_user()),
      hidden(div(id="username",radioButtons("userrev","Username?", choiceNames=list("rdmorin","HoumanLM","ninaliuta","callumcbrown","k.dreval","cmattsson","gillissierra","mannycruz","Elritch","hayashaalan","kcoyle","obigriffith"),choiceValues=list("rdmorin","HoumanLM","ninaliuta","callumcbrown","k.dreval","cmattsson","gillissierra","mannycruz","Elritch","hayashaalan","kcoyle","obigriffith"),inline = TRUE))),
      actionButton("submit", "Submit your rating!"),
      imageOutput("photo")
    )
  )
)
server <- function(input, output, session) {
  # Whenever a field is filled, aggregate all form data
  formData <- reactive({
    data <- sapply(fields, function(x) input[[x]])
    data
  })
  
  #radio button controlling how to subset the data
  subset = reactive({
    if(input$status=="reviewed"){
      show("username")
      oldrev = filter(review_results,input$userrev == review_results$user) %>% pull(mutation)
      newrev <- filter(review_results, (!review_results$mutation %in% oldrev)) %>% distinct(mutation)
      dplyr::filter(options_df,basename %in% newrev$mutation)
    }else{
      hide("username")
      dplyr::filter(options_df,!basename %in% review_results$mutation)
    }
  })
  
  observeEvent(subset(), {
    choices <- unique(subset()$Gene)
    updateSelectInput(inputId = "gene", choices = choices) 
  })
  
  g = reactive({
    subset() %>% dplyr::filter(Gene==input$gene)
  })
  observeEvent(g(), {
    choices <- unique(g()$basename)
    updateSelectInput(inputId = "mutation", choices = choices) 
  })
  
  
  
  # randomly pick a gene for the user
  observeEvent(input$random, {
    subset_df = subset()
    #limit the choices to the genes with at least one qualifying variant available
    this_gene = sample(unique(subset_df$Gene),1)
    updateSelectInput(inputId = "gene", choices = this_gene)
    
  })
  # When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    saved = save_data(formData())
    review_results <<- saved
    View(review_results)
    output$leaderboard <- renderTable({
      saved %>% group_by(user) %>% 
        tally() %>% arrange(desc(n)) %>%
        rename(c("Reviewer"="user","Submissions"="n"))
    })
    if(input$status=="unreviewed"){
      choices = g() %>%
        dplyr::filter(!basename %in% saved$mutation) %>%
        pull(basename) %>%
        unique()
    }else{
      oldrev = filter(review_results,input$userrev == review_results$user) %>% pull(mutation)
      newrev = filter(review_results, (!review_results$mutation %in% oldrev)) %>% distinct(mutation)
      choices = g() %>%
        dplyr::filter(basename %in% newrev$mutation) %>%
        pull(basename)
    }
    updateSelectInput(inputId = "mutation", choices = choices) 
    updateRadioButtons(inputId="rb", 
                       label="Mutation quality:",
                       choiceNames = list(
                         "0: Zero support for the variant in the reads",
                         "1: Minimal support and/or severe confounders",
                         "2: Low support (2-3 molecules) or other confounders",
                         "3: Modest support but some uncertainty or a confounder",
                         "4: Good support",
                         "5: Excellent support, no ambiguity"
                       ),
                       choiceValues = list("0", "1", "2","3","4","5")
    )
    updateCheckboxGroupInput(inputId="tag",
                             label="Select any tags that apply (Optional)",
                             choiceNames = c("Adjacent indel: attributable to misalignment caused by a nearby insertion or deletion",
                                             "Directional: Variant is only/mostly found on reads in the same direction",
                                             "Multiple Variants: Variant appears multi-allelic, having read support for two or more non-reference alleles",
                                             "Low Variant Frequency: Low VAF in tumour",
                                             "End of reads: Variant is only seen within ~30 bp of end of all variant-supporting reads",
                                             "High Discrepancy Region (supported by reads that have other recurrent mismatches across the track)",
                                             "Multiple Mismatches: Variant is supported by reads that have other mismatched base pairs",
                                             "Low Mapping quality: Variant is mostly supported by reads that have low mapping quality",
                                             "Short Inserts Only: Variant is exclusively found on short fragments such that sequencing from each end results in overlapping reads",
                                             "Low Count Tumor: Variant has inadequate coverage in the tumour",
                                             "Same Start End: Variant is only observed in reads that start and stop at the same positions",
                                             "No Count Normal: No coverage in the matched normal or data from matched normal is unavailable",
                                             "Mononucleotide repeat: adjacent to a region in the genome that has a single-nucleotide repeat (e.g., AAAAAA…)",
                                             "Dinucleotide repeat: adjacent to a region in the genome that has two alternating nucleotides (e.g. TGTG...)",
                                             "Tandem repeat: adjacent to a region in the genome that has three or more alternating nucleotides (e.g., GTGGTGGTG…)",
                                             "Ambiguous other: surrounded by inconclusive genomic features that cannot be explained by other tags"),
                             choiceValues = c("AI","D","MV","LVF","E","HDR","MM","LM","SIO","LCT","SSE",
                                              "NCN","MN","DN","TR",
                                              "AO"),
                             selected = "NCN",inline = F)
  })

  output$todo <- renderText({
    req(input$mutation)
    ref = options_df %>% filter(basename==input$mutation) %>% pull(Ref)
    alt = options_df %>% filter(basename==input$mutation) %>% pull(Alt)
    if(alt == "snapshot.png"){
      alt="-"
    }
    paste0("variant: ",ref,">",alt)
    
  })
  
  output$photo <- renderImage({
    req(input$mutation)
    this_row = g() %>% 
      filter(basename==input$mutation)
    full_path = pull(this_row,full)
    list(
      src = paste0("./", full_path),
      contentType = "image/png"
      #width = w,
      #height = 800
    )
  }, deleteFile = FALSE)
}
shinyApp(ui, server)
