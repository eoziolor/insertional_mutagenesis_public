# Writing function to parse cigar strings
cigar.parse<-function(obj){
  dat<-data.frame(tot.map=c(),
                  tot.idm=c())
  for(i in 1:length(obj)[1]){
    # Finding all matched reads
    eqs<-as.numeric(gsub("M","",str_extract_all(obj[i],"[0-9]+M")[[1]]))
    ins<-as.numeric(gsub("S","",str_extract_all(obj[i],"[0-9]+S")[[1]]))
    dels<-as.numeric(gsub("H","",str_extract_all(obj[i],"[0-9]+H")[[1]]))
    temp.dat<-data.frame(tot.map = sum(eqs, na.rm = TRUE),
                         tot.idm = sum(c(ins,dels),na.rm = TRUE))
    dat<-rbind(dat,temp.dat)
    if(i/1000 == round(i/1000)){
      print(i)
    }
  }
  return(dat)
}