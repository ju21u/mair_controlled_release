load_files = function(current.research.flight, res){
      #summary.current.RF <<- read.csv(file.path("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/", current.research.flight, "/5x1/all.segs.csv"), sep = " ")
      if (current.research.flight == "RF04"){
        summary.current.RF <<- read.csv(file.path("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/", current.research.flight, "/5x1/all.segs.csv"), sep = " ")
        
        current.date <<- "2021-07-30"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20211028_load_data_RF04.R"
        if (res == "1x1"){
          current.path = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/RF04/1x1/SEGMENT_RASTERS/"
        } else if (res == "5x1"){
          current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF04/5x1/SEGMENT_RASTERS/"
        }

      }else if (current.research.flight == "RF05"){
        summary.current.RF <<- read.csv(file.path("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/", current.research.flight, "/5x1/all.segs.csv"), sep = " ")
        
        current.date <<- "2021-08-03"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20211101_load_data_temp.R"
        summary.current.RF$seg <<- seq(length(summary.current.RF$start))
        summary.current.RF <<- summary.current.RF[rowSums(is.na(summary.current.RF)) == 0, ]   
    
        if (res == "1x1"){
          current.path =  "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/RF05/1x1/SEGMENT_RASTERS" 
          #"/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF05/1x1/SEGMENT_RASTERS"
        } else if (res == "5x1"){
          current.path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF05/5x1/SEGMENT_RASTERS"
        }
      } else if (current.research.flight == "RF06"){
        current.date<<-  "2021-08-06"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/load_data_RF06.R"
        summary.current.RF <<- read.csv(file.path("/n/holyscratch01/wofsy_lab/chulakadabba/RF06/", "all.segs.csv"))
        current.path = "/n/holyscratch01/wofsy_lab/chulakadabba/RF67"
        if (resolution== "5x1"){
          current.path =  "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/Level3_2023/MethaneAIR_2021/RF06/CO2_Proxy_5x1_Sax10/Level3/"
        } 
        
        
        # current.date<<-  "2021-08-06"
        # source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20230405_load_data_RF06.R"
        # #  source_file = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF06/20220213_load_data_RF06.R"
        # if (res == "1x1"){
        #   current.path =   "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/MiVida/RF06/"
        # } 
        # 
        # if (res == "5x1"){
        #   current.path =   "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/Level3_2023/MethaneAIR_2021/RF06/CO2_Proxy_5x1_Sax10/Level3/"
        #   } 
      }else if (current.research.flight == "RF07"){
        
        current.date <<- "2021-08-09"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/load_data_RF07.R"
        summary.current.RF <<- read.csv(file.path("/n/holyscratch01/wofsy_lab/chulakadabba/RF07/", "all.segs.csv"))
        current.path = "/n/holyscratch01/wofsy_lab/chulakadabba/RF07"
        if (resolution== "5x1"){
          current.path =  "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/Level3_2023/MethaneAIR_2021/RF07/CO2_Proxy_5x1_Sax10/Level3/"
        } 
        
        # current.date <<- "2021-08-09"
        # source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20230405_load_data_RF07.R"
        # 
        # #source_file = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF07/20220213_load_data_RF07.R"
        # if (res == "1x1"){
        #   current.path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/MiVida/RF07/"
        # } 
        # 
        # if (res == "5x1"){
        #   current.path =   "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/Level3_2023/MethaneAIR_2021/RF07/CO2_Proxy_5x1_Sax10/Level3/"
        # } 
        
      }else if (current.research.flight == "RF08"){
        current.date <<- "2021-08-11"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/20230124_load_data_RF08.R"
          #"/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20211028_load_data_RF04.R"
          summary.current.RF <<- read.csv(file.path("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/From_Steve/all.segs/", current.research.flight, "/5x1/all.segs.csv"))
          
          if (resolution== "5x1"){
          current.path =  "/n/holyscratch01/wofsy_lab/chulakadabba/RF08/"
        } 
        
      }else if (current.research.flight == "RF01E"){
        current.date <<- "2022-10-25"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/load_data_RF01E.R"
        summary.current.RF <<- read.csv(file.path("/n/holyscratch01/wofsy_lab/chulakadabba/MAIRE/", "RF01E.all.segs.csv"))
        current.path = "/n/holyscratch01/wofsy_lab/chulakadabba/MAIRE/OUTPUT_RASTERS/"
        
        if (resolution== "5x1"){
          current.path =  "/n/holylfs04/LABS/wofsy_lab/Lab/msargent/MAIR_E/L3/5x1_20m/CONTROLLED_RELEASES/OUTPUT_RASTERS"
        } 
        
      }else if (current.research.flight == "RF03E"){
        current.date <<- "2022-10-29"
        source_file = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/load_data_RF03E.R"
        summary.current.RF <<- read.csv(file.path("/n/holyscratch01/wofsy_lab/chulakadabba/MAIRE/", "RF03E.all.segs.csv"))
        current.path = "/n/holyscratch01/wofsy_lab/chulakadabba/MAIRE/OUTPUT_RASTERS/"
        if (resolution== "5x1"){
          current.path =   "/n/holylfs04/LABS/wofsy_lab/Lab/msargent/MAIR_E/L3/5x1_20m/CONTROLLED_RELEASES/OUTPUT_RASTERS"
        } 
        
      }
  return(list(source_file, current.path))
}