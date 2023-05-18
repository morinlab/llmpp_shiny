library(shiny)
library(tidyverse)


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
        choiceNames = c("Adjacent indel","Ambiguous other","Directional",
                        "Multiple Variants","Mononucleotide repeat","Dinucleotide repeat","Tandem repeat","Low Variant Frequency",
                        "End of reads","High Discrepancy Region","Multiple Mismatches",
                        "Low Mapping quality","Short Inserts Only",
                        "Low Count Tumor","Same Start End",
                        "Low Count Normal","No Count Normal","Tumor in Normal"),
        choiceValues = c("AI","AO","D","MV","MN","DN","TR","LVF",
                         "E","HDR","MM","LM","SIO","LCT","SSE",
                         "LCN","NCN","TN"),
        selected = NULL
        ,inline = T
      ),
    tableOutput("leaderboard")
    ),
  mainPanel(
    textInput("comment","Enter a comment (optional)"),
    
    textInput("user", "Please change this to your GitHub user ID:",value = get_user()),
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
      dplyr::filter(options_df,basename %in% review_results$mutation)
    }else{
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
    output$leaderboard <- renderTable({
      saved %>% group_by(user) %>% 
        tally() %>% arrange(desc(n)) %>%
        rename(c("Reviewer"="user","Submissions"="n"))
    })
    choices = g() %>% 
      dplyr::filter(!basename %in% saved$mutation) %>%
      pull(basename) %>%
      unique()
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
                             choiceNames = c("Adjacent indel","Ambiguous other","Directional",
                                             "Multiple Variants","Mononucleotide repeat","Dinucleotide repeat",
                                             "Tandem repeat","Low Variant Frequency",
                                             "End of reads","High Discrepancy Region",
                                             "Multiple Mismatches",
                                             "Low Mapping quality","Short Inserts Only",
                                             "Low Count Tumor","Same Start End",
                                             "Low Count Normal","No Count Normal","Tumor in Normal"),
                             choiceValues = c("AI","AO","D","MV","MN","DN","TR","LVF",
                                              "E","HDR","MM","LM","SIO","LCT","SSE",
                                              "LCN","NCN","TN"),
                             selected = NULL
                             ,inline = T)
  })
  output$todo <- renderText({
    numleft = subset() %>% nrow()
    paste("Mutations unreviewed:",numleft)
    
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
