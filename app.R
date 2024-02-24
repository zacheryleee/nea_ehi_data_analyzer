
# #For running R Shiny
library(shiny)
library(shinythemes)
library(data.table)
library(ggplot2)
library(rsconnect)

#
#
# # For the modelling (run it twice due to some issue with dependencies among the packages)
library('ISLR')
library('mixtools')
library(stats4)
library('devtools')
library(cutoff)
library('bbmle')
library("usethis")
library("dplyr")
#
# # For the Smart Gravitrap
library("lubridate")
library(tibble)
library(dplyr)
library(openxlsx)

# User Interface pages 
# about page
about_page <- tabPanel(
  title = "About", 
  titlePanel("About"),
  br(),
  h1(HTML("<b>Pupae Size Anaylsis</b>"), 
     style="text-align:left"),
  "The purpose of this web application is to accurately split a bimodal distribution to 2 normal distributions with a confidence interval of 99.9% and Type1 error of 0.1% through computational modelling.",
  br(),
  br(),
  h1(HTML("<b>Smart Gravitrap Clean-up</b>"), 
     style="text-align:left"),
  "The purpose of this web application is to clean up the raw file from the Smart Gravitrap for further anaylsis",
  br(),
  br(),
  br(),
  "Solely for use at Environmental Health Institute (EHI), NEA @ AMK Techplace II",
  br(),
  "Created by Zachery Lee Wei Quan with R Shiny",
  br(),
  "2022 August, Version I",
  br(),
  "All rights reserved"
)



# main page - Pupae Analyser
main_page1 <- tabPanel(
  title = "Pupae Size Anaylsis",
  titlePanel("Analysis"),
  sidebarLayout(
    sidebarPanel(
      title ="Inputs",
      fileInput("csv_input1", "Select Group File(s) to Import", accept = ".csv", multiple=T),
      br(),
      actionButton("run_button1", "Run Analysis", icon = icon("play"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title= "Plot",
          plotOutput("plot_1"),
          br(),
          br(),
          tableOutput("text_1")
        ),
        tabPanel(
          title = "Data File",
          fluidRow(
            column(width =12, tableOutput("combined_summary_table1"))
          )
        )
      )
    )
  )
)


# main page - Smart Gravitrap Analyzer
main_page <- tabPanel(
  title = "Smart Gravitrap Clean-up",
  titlePanel("Analysis"),
  sidebarLayout(
    sidebarPanel(
      title ="Inputs",
      fileInput("csv_input", "Select RAW Smart Gravitrap File to Import", accept = c(".csv", "text/comma-separated-values", " text/csv"), multiple=F),
      actionButton("run_button", "Run Clean-up", icon = icon("play")),
      downloadButton("download", "Download Here")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title= "Data File",
          fluidRow(
            column(width=12,dataTableOutput("combined_summary_table")),
          )
        )
      )
    )
  )
)






# user interface creation 
ui <- navbarPage(
  title = " EHI Data Anaylser",
  theme = shinytheme("sandstone"),
  main_page1,
  main_page,
  about_page
)






##############################################################
#Function - For Pupae Anaylzer 
create_combined_table1 <- function(data_input1){
  data_input1 <- unlist(data_input1, use.names=F)
  data_input1 <- data_input1[-which(is.na(data_input1))]
  length(data_input1) <- prod(dim(matrix(data_input1, ncol=12)))
  ddd <- matrix(data_input1, ncol=12 )
  colnames(ddd) <- c(paste("#",seq(1:12),sep=""))
  return(ddd)
} 



plot_mix_comps <- function(x, mu, sigma,lam){
  lam * dnorm(x, mu, sigma)
}



draw_plot_1 <- function(data_input1){
  data_input1 <- unlist(data_input1, use.names=F)
  data_input1<- data_input1[-which(is.na(data_input1))]
  
  set.seed(1234)
  mixmdl = normalmixEM(data_input1)
  mixmodel <- em(data_input1, "normal", "normal")
  confint(mixmodel,level=0.99)
  cut_off <- cutoff(mixmodel, distr =1, level=0.999, type1 = 0.001)
  
  
  data.frame(x= mixmdl$x) %>%
    ggplot() +
    geom_histogram(aes(x, ..density..), binwidth = 0.01, colour = "black", 
                   fill = "grey") + scale_x_continuous(breaks = seq(0.80, 1.45, by=0.025))  + coord_cartesian(xlim = c(0.80, 1.45))+ 
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                  aes(color = "Male"), lwd = 1, lty=1, show.legend=T, inherit.aes =F) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                  aes(color ="Female"), lwd = 1, lty=1, show.legend=T, inherit.aes=F) + 
    scale_colour_manual("Gender", values = c("pink3", "blue"), labels = c("Female", "Male")) +
    ylab("Density") + xlab("Length of Pupae (mm)") + 
    geom_vline(xintercept = cut_off[1], lty=2, lwd=1.0, col="darkblue") + 
    annotate(x=cut_off[1], y= +Inf, label = paste("Estimate Cutoff = ",round(cut_off[1],digit=3),"mm",sep=""), vjust =3, geom = "label", colour = "red", size =4.5, fill = "white" )+
    annotate(x=cut_off[1], y= +Inf, label = paste("Male Percentile = ",signif(pnorm((cut_off[1] - mixmdl$mu[1])/mixmdl$sigma[1])*100, digit=2),"%",sep=""), vjust =4.5, geom = "label", colour = "darkblue", size =4, fill ="white")+
    theme(axis.text.x =element_text(size=6.5), axis.title=element_text(size=10,face="bold.italic")) +
    theme(legend.title = element_text(face="bold.italic", size =10), legend.position = c(.95,.75))
}

write_text_1 <- function(data_input1){
  data_input1 <- unlist(data_input1, use.names=F)
  data_input1 <- data_input1[-which(is.na(data_input1))]
  
  set.seed(1234)
  mixmdl = normalmixEM(data_input1)
  mixmodel <- em(data_input1, "normal", "normal")
  confint(mixmodel,level=0.99)
  cut_off <- cutoff(mixmodel, distr =1, level=0.999, type1 = 0.001)
  
  testing <-as.data.frame(mixmodel$param)
  testing[5,] <- pnorm((cut_off[1] - mixmdl$mu[1])/mixmdl$sigma[1])*100
  testing <- t(testing)
  colnames(testing) <- c("Male Pupae Mean (mm)", "Male Pupae Standard Deviation (mm)","Female Pupae Mean (mm)", "Female Pupae Standard Deviation (mm)", "Male Pupae Percentile Cutoff (%)")
  
  return(testing)
} 






#Function - For Smart Gravitrap Analyzer 
#between function
between <- function(x,lower,upper,incbounds=T)
{
  if(T) x>=lower & x<=upper
  else (F) 
} %>% "%between%" <- function(x,y) between(x,y[1],y[2],incbounds=T)


# function to see the workbook as a dataframe
create_combined_table <- function(x){
  
  # for making of the data file for viewing, unable to download yet 
  testingwb1 <- createWorkbook()
  addWorksheet(testingwb1, "Smart Gravitrap1")
  headerStyle1 <- createStyle(fontSize = 13, fontColour = "black",
                              border = "TopBottomLeftRight", borderColour = "black",
                              textDecoration = "bold", wrapText = F)
  
  writeData(testingwb1, "Smart Gravitrap1", "Date", startCol=1, startRow = 1) +
    writeData(testingwb1, "Smart Gravitrap1", "E-Week", startCol=2, startRow = 1) +
    writeData(testingwb1, "Smart Gravitrap1", "Time-interval", startCol=3, startRow = 1) +
    writeData(testingwb1, "Smart Gravitrap1", "Start Time", startCol=4, startRow = 1) +
    writeData(testingwb1, "Smart Gravitrap1", "End Time ", startCol=5, startRow = 1) +
    writeData(testingwb1, "Smart Gravitrap1", "Label", startCol=6, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", "Time", startCol=7, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", "Remark", startCol=8, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", " M.Ae.aegypti", startCol=9, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", "F.Ae.aegypti", startCol=10, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", " M.Ae.albopictus", startCol=11, startRow= 1) +
    writeData(testingwb1, "Smart Gravitrap1", " F.Ae.albopictus", startCol=12, startRow= 1) 
  
  addStyle(testingwb1, "Smart Gravitrap1", headerStyle1, rows=1, cols=seq(1:12), gridExpand =T)
  
  testing_curr_row <- 2
  for (viii in seq_along(x)){
    writeData(testingwb1, "Smart Gravitrap1", names(x)[viii], startCol=1, startRow = testing_curr_row)
    writeData(testingwb1, "Smart Gravitrap1", strftime(names(x)[viii], format ="%W"), 
              startCol=2, startRow = testing_curr_row)
    for (ix in names(x[[viii]])){
      writeData(testingwb1, "Smart Gravitrap1", ix, startCol=3, startRow = testing_curr_row)
      writeData(testingwb1, "Smart Gravitrap1", gsub(" ", "", strsplit(ix, split = "-")[[1]][1]), startCol=4, startRow = testing_curr_row)
      writeData(testingwb1, "Smart Gravitrap1", gsub(" ", "", strsplit(ix, split = "-")[[1]][2]), startCol=5, startRow = testing_curr_row)
      writeData(testingwb1, "Smart Gravitrap1", x[[viii]][[which(names(x[[viii]]) == ix)]][,c(-2,-3)], 
                startCol=6, startRow= testing_curr_row, colNames=F)
      
      writeData(testingwb1, "Smart Gravitrap1", rep(names(x)[viii], nrow(x[[viii]][[which(names(x[[viii]]) == ix)]])), startCol=1, startRow = testing_curr_row )
      writeData(testingwb1, "Smart Gravitrap1", rep(strftime(names(x)[viii], format ="%W"), nrow(x[[viii]][[which(names(x[[viii]]) == ix)]])), startCol=2, startRow = testing_curr_row)
      writeData(testingwb1, "Smart Gravitrap1", rep(ix, nrow(x[[viii]][[which(names(x[[viii]]) == ix)]])), startCol=3, startRow = testing_curr_row)
      writeData(testingwb1, "Smart Gravitrap1", rep(gsub(" ", "", strsplit(ix, split = "-")[[1]][1]), nrow(x[[viii]][[which(names(x[[viii]]) == ix)]])), startCol=4, startRow = testing_curr_row )
      writeData(testingwb1, "Smart Gravitrap1", rep(gsub(" ", "", strsplit(ix, split = "-")[[1]][2]), nrow(x[[viii]][[which(names(x[[viii]]) == ix)]])), startCol=5, startRow = testing_curr_row )
      
      if (viii %% 2 == 0){
        addStyle(testingwb1, "Smart Gravitrap1", createStyle(fgFill= "lightcyan3", border = "TopBottomLeftRight"), rows=c(testing_curr_row: (testing_curr_row+nrow(x[[viii]][[which(names(x[[viii]]) == ix)]]))), cols=seq(1:12), gridExpand =T)
      } else {
        addStyle(testingwb1, "Smart Gravitrap1", createStyle(fgFill= "white", border = "TopBottomLeftRight"), rows=c(testing_curr_row: (testing_curr_row+nrow(x[[viii]][[which(names(x[[viii]]) == ix)]]))), cols=seq(1:12), gridExpand =T)
      }
      testing_curr_row <- testing_curr_row + nrow(x[[viii]][[which(names(x[[viii]]) == ix)]]) 
    }
  }
  testing <- readWorkbook(testingwb1, sheet=1)
  return(testing)
}



# function to convert to clean the data_input
clean_data <- function(y){
  #Important column to take note 
  x<- c("Label","Submitted Date", "Submitted Time", "Remark", "M.Ae.aegypti", "F.Ae.aegypti", "M.Ae.albopictus", "F.Ae.albopictus")
  
  index <- which(colnames(y) %in% x)
  
  trap1 <- y[,..index]
  
  
  # Converting the time to a six numerical figure for easy conversion later
  temp566 <- c()
  for (i in trap1[,3]){
    for (ii in i){
      while (nchar(ii) != 6){
        ii <- paste0("0",ii) 
      } 
      temp566 <- append(temp566, format(strptime(ii, format = "%H%M%S"),"%H:%M:%S"))
    }
    trap1[,3] <- temp566
  }
  
  
  # changing the number in date column to date format and clean up any unusual dates
  trap1[,2] <- sapply(trap1[,2], function(x) as.character(as.Date.character(x, format= "%Y%m%d")))
  #trap1 <- trap1[trap1$Submitted.Date >= "2021-01-01",]
  
  ### Cleaning the label#####
  trap1$Label <- sapply(trap1$Label, function(x) paste("T0000", x, sep=""))
  
  #### Start of iterative codes
  
  trap1 <-trap1[order(trap1[,2]),]
  
  d <- table(trap1[,2])
  list1<- list()
  for (dd in row.names(d)) {
    list1[dd] <- list(NULL)
  }
  
  # sorting the dates together for  analysis
  for (ii in names(list1)){
    tempindex <- which(trap1[,2] == ii)
    list1[[ii]] <- trap1[c(tempindex),]
  }
  for (iii in 1:length(list1)){
    Time_interval <- NA
    list1[[iii]] <- add_column(list1[[iii]], Time_interval, .after= 2)
    
    list1[[iii]][,3][list1[[iii]][,4] %between% c("00:00:00", "00:59:59")] <- "12am - 1am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("01:00:00", "01:59:59")] <- "1am - 2am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("02:00:00", "02:59:59")] <- "2am - 3am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("03:00:00", "03:59:59")] <- "3am - 4am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("04:00:00", "04:59:59")] <- "4am - 5am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("05:00:00", "05:59:59")] <- "5am - 6am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("06:00:00", "06:59:59")] <- "6am - 7am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("07:00:00", "07:59:59")] <- "7am - 8am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("08:00:00", "08:59:59")] <- "8am - 9am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("09:00:00", "09:59:59")] <- "9am - 10am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("10:00:00", "10:59:59")] <- "10am - 11am"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("11:00:00", "11:59:59")] <- "11am - 12pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("12:00:00", "12:59:59")] <- "12pm - 1pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("13:00:00", "13:59:59")] <- "1pm - 2pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("14:00:00", "14:59:59")] <- "2pm - 3pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("15:00:00", "15:59:59")] <- "3pm - 4pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("16:00:00", "16:59:59")] <- "4pm - 5pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("17:00:00", "17:59:59")] <- "5pm - 6pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("18:00:00", "18:59:59")] <- "6pm - 7pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("19:00:00", "19:59:59")] <- "7pm - 8pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("20:00:00", "20:59:59")] <- "8pm - 9pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("21:00:00", "21:59:59")] <- "9pm - 10pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("22:00:00", "22:59:59")] <- "10pm - 11pm"
    list1[[iii]][,3][list1[[iii]][,4] %between% c("23:00:00", "23:59:59")] <- "11pm - 12am"
    
  }
  
  ## to split the dataframe base on the time interval 
  for (iv in 1:length(list1)){
    list1[[iv]] <- split(list1[[iv]], factor(list1[[iv]]$Time_interval, levels=unique(list1[[iv]]$Time_interval)))
  }
  
  return(create_combined_table(list1))
}


# Server 
server <- function(input, output){
  options(shiny.maxRequest = 20*1024^2)
  
  # For Pupae Anaylzer 
  data_input1 <- reactive ({
    req(input$csv_input1)
    rbindlist(lapply(input$csv_input1$datapath, fread), use.names =F, fill= F)
  })
  
  # for the statistic table 
  combined_summary_table1 <- eventReactive(input$run_button1,{
    create_combined_table1 (data_input1())
  })
  output$combined_summary_table1 <- renderTable(combined_summary_table1(), bordered =T, width = "100%")
  
  #for the plot table
  plot_1 <- eventReactive(input$run_button1,{
    draw_plot_1(data_input1())
  })
  output$plot_1 <- renderPlot(plot_1(), width = "auto", height ="auto")
  
  # for the text 
  text_1 <- eventReactive(input$run_button1,{
    write_text_1(data_input1())
  })
  output$text_1 <- renderTable(text_1())
  
  # For Smart Gravitrap Anaylzer 
  data_input <- reactive({
    req(input$csv_input)
    fread(input$csv_input$datapath, header=T)
  })
  
  
  # for the statstic table 
  combined_summary_table <- eventReactive(input$run_button,{
    clean_data(data_input())
  })
  output$combined_summary_table <- renderDataTable(combined_summary_table())
  
  
  # For the download button 
  output$download <- downloadHandler(
    filename = function() {
      paste(input$csv_input,"CLEANED.csv", sep="_")
    },
    content = function(file) {
      write.csv(combined_summary_table(), file, row.names=F)
    }
  )
}



shinyApp(ui = ui, server = server)

