ui <- navbarPage("methylscaper",
                 tabPanel("Single-cell",
                          navbarPage("",
                                     tabPanel("Seriation",
                                              sidebarLayout(
                                                  sidebarPanel(
                                                      fileInput("gch_seq_file",
                                                                label="GCH 
                                                                Sequence RDS 
                                                                file"),
                                                      fileInput("hcg_seq_file",
                                                                label="HCG 
                                                                Sequence RDS 
                                                                file"),
                                                      numericInput("startPos", 
                                                                   label=
                                                                   "Starting 
                                                                   position", 
                                                                   value=0, 
                                                                   min=1),
                                                      numericInput("endPos", 
                                                                   label=
                                                                   "Ending 
                                                                   position", 
                                                                   value=0, 
                                                                   min=1),
                                                      uiOutput
                                                      ("positionSlider"),
                                                      selectInput(
                                                                  "sc_ser_
                                                                  method", 
                                                                  label=
                                                                  "Seriation 
                                                                  Method:",
                                                                  choices=
                                                                  c("PCA", 
                                                                  "ARSA")),
                                                      selectInput("sc_refine_
                                                                  method", 
                                                                  label=
                                                                  "Refinement 
                                                                  Method:",
                                                                  choices=
                                                                  c("PCA", 
                                                                    "HC_
                                                                    average")),
                                                      radioButtons("sc_brush_
                                                                   choice", 
                                                                   label=
                                                                   "Brushing 
                                                                   for:",
                                                                   choices=c(
                                                                   "Refinement"
                                                                   , "Weight
                                                                   ing"), 
                                                                   selected=
                                                                   "Weighting"
                                                                   ),
                                                      actionButton(
                                                                   "sc_force_
                                                                   reverse", 
                                                                   label="Force
                                                                    Reverse"),
                                                      verbatimTextOutput
                                                      ("sc_info")
                                                ),

                                                mainPanel(
                                                  fluidRow(column(width=8,
                                                                  plotOutput
                                                                  (outputId=
                                                                   "sc_seqPlot"
                                                                    , brush=
                                                                    "sc_plot_
                                                                    brush",  
                                                                    width=
                                                                    "100%")),
                                                           column(width=2, 
                                                                  align='left',
                                                                  selectInput
                                                                  ("sc_plot_
                                                                   filetype", 
                                                                   label="File 
                                                                   type", 
                                                                   choices=
                                                                   c("PNG", 
                                                                     "SVG", 
                                                                     "PDF")),
                                                            downloadButton
                                                            ("sc_plot_down", 
                                                            label="Download
                                                            the plot"),
                                                            downloadButton
                                                            ("sc_log_down", 
                                                            label = "Download 
                                                            changes log")
                                                            )

                                                  )
                                                )
                                              )),
                                     tabPanel("Summary Statistics",
                                              radioButtons("sc_
                                                           proportion_choice", 
                                                           label="Proportion 
                                                           of:", choices=
                                                           c("Yellow", "Red"), 
                                                           selected="Yellow"),
                                              splitLayout(cellWidths=c("50%",
                                                                       "50%"),
                                                  plotOutput(outputId="sc_
                                                             proportion_color_
                                                             histogram"),
                                                  plotOutput(outputId="sc_
                                                             percent_C")),
                                              splitLayout(cellWidths=c("50%", 
                                                                       "50%"),
                                                downloadButton("sc_proportion_
                                                               hist_download", 
                                                               label="Download 
                                                               histogram"),
                                                downloadButton("sc_percentC_
                                                               plot_download", 
                                                               label="Download 
                                                               plot")),
                                              splitLayout(cellWidths=c("50%", 
                                                                       "50%"),
                                                downloadButton("sc_proportion_
                                                               data_download", 
                                                               label="Download 
                                                               proportion 
                                                               data"),
                                                downloadButton("sc_percentC_
                                                               data_download", 
                                                               label="Download 
                                                               percentage data"
                                                               ))))),
                 tabPanel("Single-molecule",
                          navbarPage("",
                                     tabPanel("Preprocessing",
                                              fileInput("fasta.file", 
                                                        label="FASTA File"),
                                              fileInput("ref.file", 
                                                        label="Reference 
                                                        File"),
                                              textInput("gch.file.name", 
                                                        label="GCH File Name"),
                                              textInput("hcg.file.name", 
                                                        label="HCG File Name"),
                                              textInput("processing.log.name", 
                                                        label="Processing Log 
                                                        File Name"),
                                              actionButton("run.align", 
                                                           label="Run")),
                                     tabPanel( "Seriation",
                                               sidebarLayout(
                                                 sidebarPanel(
                                                   fileInput("sm_gch_file", 
                                                             label="GCH Data 
                                                             file input"),
                                                   fileInput("sm_hcg_file", 
                                                             label="HCG Data 
                                                             file input"),
                                                   selectInput("sm_ser_method", 
                                                               label="Seriation
                                                                Method:",
                                                               choices=c("PCA", 
                                                                         "ARSA"
                                                                         )),
                                                   selectInput("sm_refine_
                                                               method", 
                                                               label=
                                                               "Refinement 
                                                               Method:",
                                                               choices=c("PCA", 
                                                                         "HC_
                                                                         avera
                                                                         ge")),
                                                   radioButtons("sm_brush_
                                                                choice", 
                                                                label="Brushing 
                                                                for:",
                                                                choices=
                                                                c("Refinement", 
                                                                  "Weighting"), 
                                                                  selected=
                                                                  "Weighting"),
                                                   actionButton("sm_force_
                                                                reverse", 
                                                                label="Force 
                                                                 Reverse"),
                                                   verbatimTextOutput
                                                   ("sm_info")
                                                 ),

                                                 mainPanel(
                                                   fluidRow(column(width=8,
                                                                   plotOutput
                                                                   (outputId=
                                                                    "sm_
                                                                    seqPlot",
                                                                    brush="sm_
                                                                    plot_brush"
                                                                    ,  width=
                                                                    "100%")),
                                                            column(width=2, 
                                                                   align='left',
                                                                   selectInput
                                                                   ("sm_
                                                                    filetype",
                                                                    label="File
                                                                    type", 
                                                                    choices=
                                                                    c("PNG", 
                                                                      "SVG", 
                                                                      "PDF")),
                                                              downloadButton("sm_
                                                                          plot_
                                                                          down"
                                                                          , 
                                                                          label
                                                                          ="Dow
                                                                          nload 
                                                                          the 
                                                                          plot"
                                                                          ),
                                                              downloadButton
                                                              ("sm_log_down", 
                                                                label="Downl
                                                                oad changes
                                                                log")
                                                              )

                                                   )
                                                 )
                                               )),
                                     tabPanel("Summary Statistics",
                                              radioButtons("sm_proportion_
                                                           choice", 
                                                           label="Proportion 
                                                           of:", choices=
                                                          c("Yellow", "Red"),
                                                          selected="Yellow"),
                                              splitLayout(cellWidths=
                                                          c("50%", "50%"),
                                                  plotOutput(outputId=
                                                            "sm_proportion_
                                                            color_histogram"),
                                                  plotOutput(outputId=
                                                             "sm_percent_C")),
                                              splitLayout(cellWidths=c("50%", 
                                                                       "50%"),
                                                downloadButton("sm_proportion_
                                                               hist_download", 
                                                               label="Download 
                                                               histogram"),
                                                downloadButton("sm_percentC_
                                                               plot_download", 
                                                               label="Download 
                                                               plot")),
                                              splitLayout(cellWidths=c("50%",
                                                                       "50%"),
                                                downloadButton("sm_proportion_
                                                               data_download", 
                                                               label="Download 
                                                               proportion 
                                                               data"),
                                                downloadButton("sm_percentC_
                                                               data_download", 
                                                               label="Download 
                                                               percentage data"
                                                               ))))))
