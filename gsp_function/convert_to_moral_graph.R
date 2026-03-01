convert_to_moral_graph<-function(A){
  for(i in 1:nrow(A)){
    A_0<-A
    if(sum(A_0[,i])>1){
      list1<-which(A_0[,i]==1)
      comb<-combn(list1,2)
      for(j in 1:ncol(comb)){
        A[comb[1,j],comb[2,j]]<-1
        A[comb[2,j],comb[1,j]]<-1
      }
    }
  }
  A<- t(A)+A
  A[A>1]<-1
  return(A)
}





