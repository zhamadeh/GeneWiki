library(shiny)
library(readxl)
library(dplyr)
library(DT)


myeloid=read.table("https://github.com/zhamadeh/GeneWiki/raw/main/Shiny/myeloid.txt",header=T,sep="\t",fill=T)
myeloid = subset(myeloid, !is.na(Gene.name) & Gene.name != "")
myeloid <- myeloid %>%
  arrange(desc(Consensus.Score))
colnames(myeloid)
#myeloid = myeloid[1:189,]

genes= sort(filter(myeloid,Description != "Structural aberration")$Gene.name)

abs = filter(myeloid,Description == "Structural aberration")$Gene.name

myeloid = dplyr::select(myeloid,c("Gene.name","Tier" ,"Disease","Prognosis","WHO.classification"                 
                                  ,"ICC.classification" ,"FDA.approved.drugs..OncoKB."  , "Most.common.fusion" ,"Partner.genes..variants","GSC.Myeloid.Panel" ,"UHN.Myeloid.Panel","ICC.WHO."
))
colnames(myeloid)=c("Gene","Tier" ,"Disease","Interpretation","WHO Classification"                 
                    ,"ICC Classification" ,"FDA Drugs (OncoKB)"  , "Common Fusion Partners" ,"Rare Fusion Partners","GSC Myeloid Panel" ,"UHN Myeloid Panel","ICC-WHO Publication")
#myeloid = myeloid[1:189,]



ui <- fluidPage(
  titlePanel(title = "OGM Gene Wiki"),
  tags$head(tags$style(HTML("
    .custom-sidebar {
      width: 400px; /* Set the fixed width for the sidebar */
    }
  "))),                          
  sidebarLayout(
    div( style = "width: 700px;", 
         sidebarPanel(
           div(style = "max-height: 300px; overflow-y: auto;",
               checkboxGroupInput("gene_filter", "Select Gene:",
                                  choices = c("All Genes", genes),
                                  selected = "All Genes")),
           div(style = "max-height: 300px; overflow-y: auto; margin-top: 40px;", 
               checkboxGroupInput("sv_filter", "Select Aberration:",
                                  choices = c("All Aberrations", abs),
                                  selected = "All Aberrations"))
         ),
         
    ),
    mainPanel(
      dataTableOutput("filtered_table")
    )
  )
)

server <- function(input, output) {
  
  # Create a reactive filtered dataset based on selected gene name
  filtered_data <- reactive({
    if (!("All Genes" %in% input$gene_filter) & !("All Aberrations" %in% input$sv_filter)) {
      dplyr::filter(myeloid, Gene %in% c(input$gene_filter, input$sv_filter))
    } else if (!("All Genes" %in% input$gene_filter) & ("All Aberrations" %in% input$sv_filter)) {
      dplyr::filter(myeloid, Gene %in% c(input$gene_filter,abs))
    } else if (("All Genes" %in% input$gene_filter) & !("All Aberrations" %in% input$sv_filter)) {
      dplyr::filter(myeloid, Gene %in% c(input$sv_filter,genes))
    } else {
      myeloid  # Show the full dataset if no filter is selected
    }
  })
  
  output$filtered_table <- DT::renderDT(server = TRUE,{
    datatable(
      filtered_data(),
      options = list(processing = TRUE,
                     autoWidth = TRUE,
                     scrollX = TRUE,
                     columnDefs = list(list(width = '10px', targets = c(1)),
                                       list(width = '10px', targets = c(2)),
                                       list(width = '100px', targets = c(3)),
                                       list(width = '400px', targets = c(4)),
                                       list(width = '200px', targets = c(5)),
                                       list(width = '200px', targets = c(6)),
                                       list(width = '10px', targets = c(7)),
                                       list(width = '10px', targets = c(8)))  # Set width for Prognosis column
      )
    )
    
  })
  
}

shinyApp(ui, server)

