
saveSet <- function(obj=NA, set="", output=1){

    if(globalConfig$anal.mode == "api"){ 
      qs:::qsave(obj, paste0(set, ".qs"));
      return(output)
    }else{
      if(set == "dataSet"){
        dataSet <<- obj;
      }else if(set == "analSet"){
        analSet <<- obj;
      }else if(set == "imgSet"){
        imgSet <<- obj;
      }else if(set == "paramSet"){
        paramSet <<- obj;
      }else if(set == "msgSet"){
        msgSet <<- obj;
      }else if(set == "cmdSet"){
        cmdSet <<- obj;
      }

        if(globalConfig$anal.mode == "web"){
            return(output);
        }else{
            return(obj);
        }
    }

}

readSet <- function(obj=NA, set=""){
    if(globalConfig$anal.mode == "api"){
      path <- "";
      if(exists('user.path')){
        path <- user.path;
      }

      if(path != ""){
        obj <- load_qs(paste0(path, set, ".qs"));
      }else{
        obj <- qs:::qread(paste0(set, ".qs"));
      }
    }
    return(obj);
}

load_qs <- function(url) qs::qdeserialize(curl::curl_fetch_memory(url)$content)

readDataset <- function(fileName=""){

    if(globalConfig$anal.mode == "api"){
      if(exists('user.path')){
        path <- user.path;
        obj <- load_qs(paste0(path, fileName));
      }else{
        obj <- qs:::qread(fileName);
      }
    }else{
       obj <- dataSets[[fileName]];
    }

    return(obj);
}