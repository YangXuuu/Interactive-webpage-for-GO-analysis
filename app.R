library(shiny)
require(tidyverse)
library(magrittr)
library(dplyr)
library(readxl)
library(ggrepel)
library(EnhancedVolcano)
library(raster)
library(stringr)
library(readr)
library(BiocManager)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

CCLE_expression <- data.table::fread("/home/yang/Project_shiny_Depmap/CCLE_expression.csv", na = "NA")
CRISPR_gene_effect <- data.table::fread("/home/yang/Project_shiny_Depmap/CRISPR_gene_effect.csv", na="NA")
colnames(CRISPR_gene_effect)[1]<-c("ID")

RNAi_database <- data.table::fread("/home/yang/Project_shiny_Depmap/D2_combined_gene_dep_scores.csv", na = "NA")
RNAi_database<- t(RNAi_database)
colnames(RNAi_database)<- RNAi_database[1,]
RNAi_database<- RNAi_database[-1,]
RNAi_database<- as.data.frame(RNAi_database)

sample_info_1_ <- read_csv("/home/yang/Project_shiny_Depmap/sample_info.csv", na = "NA")
colnames(CCLE_expression)[1] <- c("ID")

duplicated(sample_info_1_$lineage)
sample_info_1_$lineage[duplicated(sample_info_1_$lineage)]
clean1 <- sample_info_1_$lineage
clean2 <- filter(sample_info_1_,!duplicated(clean1))
clean3 <- clean2$lineage

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Differential Gene Expression Analysis"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "p_value",
                        label = "Required p value ",
                        selected = "0.01",
                        choices = c("0.01","0.05")),
            
            selectInput(inputId = "database",
                        label = "Chosen database",
                        selected = "CRISPR",
                        choices = c("CRISPR","RNAi")),
            
            textInput(inputId = "Use_gene",
                      label = "Selected gene",
                      value = "WEE1"),
            
            selectInput(inputId = "cell_line",
                      label = c(),
                      selected = "ALL",
                      choices = c("ALL",clean3)),
            
            submitButton("Submit")
                        
        ),

        # Show the generated analysis plot 
        mainPanel(
           plotOutput("volcanoPlot"),
           plotOutput("barPlot")
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    databaseInput <- reactive({
        if(input$database=="CRISPR"){
            Use_data<-CRISPR_gene_effect
            colnames(Use_data) <- gsub(".\\(.*\\)","",colnames(Use_data))
        }
        if(input$database=="RNAi"){
            Use_data<-RNAi_database
            colnames(Use_data) <- gsub(".\\(.*\\)","",colnames(Use_data))
        }
        Use_data
    })
    
    cutoffInput <- reactive({
        ## 设p value & FC cutoff阈值
        ## pCutoff = 0.01 
        if(input$p_value==0.05){
            pCutoff<-0.05
        }
        if(input$p_value==0.01){
            pCutoff<-0.01
        }
        FCcutoff =1
        c(pCutoff,FCcutoff)
    })
    
    dataInput <- reactive({
        Use_data<-databaseInput()
        Use_gene<-input$Use_gene
        ## 挑出目的基因的数据 (ID, cell line)
        select_gene<-Use_data %>%
            select(ID, Use_gene)
            colnames(select_gene)[2] <- c("READING")
        
        select_gene <- na.omit(select_gene)
    
        if(input$cell_line=="ALL"){
            select_lineage<-sample_info_1_ %>%
            select(DepMap_ID, lineage)
        }else{
            select_lineage<-sample_info_1_ %>%
                select(DepMap_ID, lineage) %>%
                filter(lineage == input$cell_line)
        }
    
    
    ## 找交集 
        vector_select_gene <- as.vector(as.matrix(select_gene[,1]))
        vector_select_lineage <- as.vector(as.matrix(select_lineage[,1]))
        intersection<- intersect(vector_select_gene,vector_select_lineage)
        select_gene <- select_gene[select_gene$ID %in% intersection]
    
    ## 确定前后5%界限 (若挑出的gene小于100个,则直接将前5个作为top&bottom)
       if(count(select_gene) <= 100){
            gene_order <- select_gene[order(select_gene$READING, decreasing = TRUE)]
            select_top <- rbind(gene_order[1],gene_order[2],gene_order[3],gene_order[4],gene_order[5])
         select_bottom <- rbind(gene_order[dim(gene_order)[1]-4],gene_order[dim(gene_order)[1]-3],gene_order[dim(gene_order)[1]-2],gene_order[dim(gene_order)[1]-1],gene_order[dim(gene_order)[1]])
       }else{
            boundary_top <- quantile(as.numeric(select_gene$READING),prob=1-95/100)
            boundary_bottom <- quantile(as.numeric(select_gene$READING),prob=1-5/100)
        ## 筛选前后5%表达
            select_top<-select_gene %>%
                dplyr::select('ID', READING) %>%
                dplyr::filter(READING < boundary_top) 
            select_bottom<-select_gene %>%
                dplyr::select('ID', READING) %>%
                dplyr::filter(READING > boundary_bottom)
        }
    
        ## 转为vector（向量）- top
        vector_top <- as.vector(as.matrix(select_top[,1]))
    
        ## 转为vector（向量）- bottom
        vector_bottom <- as.vector(as.matrix(select_bottom[,1]))
    
        ## 筛选行
        merging_top1 <- CCLE_expression[CCLE_expression$ID %in% vector_top,]
        merging_bottom1 <- CCLE_expression[CCLE_expression$ID %in% vector_bottom,]
    
       ## 把ID作为rowname，融合top&bottom数据库
        rownames(merging_top1) = merging_top1$ID
        merging_top2 <- merging_top1[,-1]
        rownames(merging_bottom1) = merging_bottom1$ID
        merging_bottom2 <- merging_bottom1[,-1]
        topbottom <- rbind(merging_top2, merging_bottom2)
    
        ## 循环执行T test
        gene_test <- as.data.frame(matrix(0,length(merging_bottom2),3))
        for(i in 1:length(merging_bottom2)){
            a <- na.omit(as.matrix(merging_top2[,..i]))
            b <- na.omit(as.matrix(merging_bottom2[,..i]))
            t <- t.test(a,b)
            t$p.value->gene_test[i,1]
            (2^mean(a)-1)/(2^mean(b)-1)->gene_test[i,2]
            colnames(merging_bottom2)[i]->gene_test[i,3]
        }
    
        ##
        na.omit(gene_test)->x
        x[x[,2]!="Inf" & x[,2]!="0",]->x
        log2(as.numeric(x[,2]))->x[,2]
        str_extract(x[,3], "\\w+")->x[,4]
        colnames(x)<-c("pvalue","log2FoldChange","gene","gene_adj")
        row.names(x)<-x[,3]
        as.numeric(x$log2FoldChange)->x$log2FoldChange
        x
    })
    
    ## 火山图
    output$volcanoPlot <- renderPlot({
        x<-dataInput()
        c<-cutoffInput()
        pCutoff<-c[1]
        fcCutoff<-c[2]
        EnhancedVolcano(x,
                    lab= x[,4],
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = pCutoff,
                    FCcutoff = fcCutoff,
                    labSize=3,
                    title="Sensitve vs Insensitve")
    })
    
    output$barPlot <- renderPlot({
        x<-dataInput()
        c<-cutoffInput()
        pCutoff<-c[1]
        FCcutoff<-c[2]
        ## 获取上调基因
        select_x<-x %>%
            dplyr::select(pvalue, log2FoldChange, gene_adj) %>%
            dplyr::filter(pvalue < pCutoff & log2FoldChange > FCcutoff)
    
        ## 获取下调基因
        select_y<-x %>%
            dplyr::select(pvalue, log2FoldChange, gene_adj) %>%
            dplyr::filter(pvalue < pCutoff & log2FoldChange < FCcutoff)
    
       #单列基因名文件
        data <-as.data.frame(select_y$gene_adj)   
        #需要character格式，然后进行ID转化
        colnames(data) <- c('V1')
        data$V1<- as.character(data$V1) 
        #将SYMBOL格式转为ENSEMBL和ENTERZID格式 
        test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
        #head(test1,2)
    
        #富集分析的背景基因集
        kk <- enrichKEGG(gene = test1$ENTREZID, minGSSize = 10,
                     organism = 'hsa', #KEGG可以用organism = 'hsa'
                     pvalueCutoff = 1,qvalueCutoff=1)
        #head(kk)
    
        ##可视化 - 条形图 - 按p从小到大排，绘制前20个Term
        barplot(kk, showCategory=20,title="EnrichmentGO_down_regulation")    
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
