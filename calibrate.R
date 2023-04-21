library(shiny)
library(tidyverse)

shiny_log <- "shiny_responses.tsv"
fields <- c("mutation","Rating","User","tag","comment")
review_results = suppressMessages(read_tsv(shiny_log,col_names=fields))

options_df = data.frame(full=dir(recursive=T,pattern=".png")) %>% 
  mutate(base=basename(full)) %>%
  filter(base %in% review_results$mutation) %>%
  mutate(basename=base) %>% 
  separate(base,into=c("Region","End","Gene","sample_id","Ref","Alt"),sep="-+") %>% 
  arrange(Gene)

seen_mutations = c()

ui <- fluidPage(
  titlePanel("See how others have rated variants"),
  sidebarLayout(
    sidebarPanel(
     actionButton("random", "Pick a random example"),
     selectInput("mutation", "Revisit a mutation", choices = unique(options_df$basename)),
     checkboxGroupInput("rb", "Rating system:",width="100%",
                        selected=NULL,
                  choices = list(
                    "0: Zero support for the variant in the reads",
                    "1: Minimal support and/or severe confounders",
                    "2: Low support (2-3 molecules) or other confounders",
                    "3: Modest support but some uncertainty or a confounder",
                    "4: Good support",
                    "5: Excellent support, no ambiguity"
                  )
     )
    ),
    mainPanel(
      tableOutput("stats"),
      imageOutput("photo")
    )
  )
)
server <- function(input, output, session) {
  
  
  observeEvent(input$view, {
    #limit the choices to the genes with at least one qualifying variant available
    base_name = input$mutation
    updateSelectInput(inputId = "mutation", choices = seen_mutations,selected = base_name) 
    output$stats <- renderTable({
      mutation_stats()
    })
    output$photo <- renderImage({
      req(input$mutation)
      
      full_path = this_image()
      #this seems to be running even if input$mutation is not set
      #causing this error:
      #'raw = FALSE' but './' is not a regular file
      list(
        src = paste0("./", full_path),
        contentType = "image/png",
        width = 1000,
        height = 800
      )
    }, deleteFile = FALSE)
  })
  
  # randomly pick a gene for the user
  observeEvent(input$random, {
    #limit the choices to the genes with at least one qualifying variant available
    mut = slice_sample(review_results,n=1) %>% pull(mutation)
    seen_mutations <<- c(mut,seen_mutations)
    print(seen_mutations)
    updateSelectInput(inputId = "mutation", choices = seen_mutations,selected = mut) 
    output$stats <- renderTable({
      mutation_stats()
    })
    output$photo <- renderImage({
      req(input$random)
      
      full_path = this_image()
      #this seems to be running even if input$mutation is not set
      #causing this error:
      #'raw = FALSE' but './' is not a regular file
      list(
        src = paste0("./", full_path),
        contentType = "image/png",
        width = 1000,
        height = 800
      )
    }, deleteFile = FALSE)
  })
  
  this_image = reactive({
    dplyr::filter(options_df,basename==input$mutation) %>% slice_sample(n=1) %>% pull(full)
  })
  mutation_stats = reactive({
    dplyr::filter(review_results,mutation==input$mutation) %>% select(-1) 
  })

}
shinyApp(ui, server)
