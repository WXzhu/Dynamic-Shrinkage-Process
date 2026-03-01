topo_order<-function(g,G,n,p,h){
  G<-unname(G)
  n_p<-colSums(G)
  col_G<-rowSums(G)
  
  n_p0<-n_p
  G0<-G
  id<-rep(0,n)
  c<-1
  if(is.dag(g)){
    while(!identical(n_p0,rep(-1,n))){
      for (i in which(n_p0==0)) {
        G0[i,]<-rep(0,p)
        G0[i,i]<--1
        id[i]<-c
      }
      c<-c+1
      n_p0<-colSums(G0)
    }
  }else{
    print("This is not a DAG.")
  }
  topo<-order(id)
  max_id<-max(id)
  
  topo_h<-topo
  level_i_nodes <- sapply(1:h, function(i) which(id == i))
  topo_h[topo_h %in% unlist(level_i_nodes)] <- 0
  
  topo_no1<-topo
  level_1_nodes <- which(id == 1)
  topo_no1[topo_no1 %in% level_1_nodes] <- 0
  
  topo_1<-topo-topo_no1
  
  topo_2<-topo
  level_2_nodes<-which(id==2)
  topo_2[!(topo_2 %in% level_2_nodes)]<-0
  
  return(list(id=id,max_id=max_id,n_p=n_p,topo=topo,topo_1=topo_1,topo_no1=topo_no1,topo_h=topo_h,topo_2=topo_2))
}
  



## b.每个节点的父节点和n_p
# n_p<-colSums(G)
# col_G<-rowSums(G)
# ## c.节点拓扑排序id
# n_p0<-n_p
# G0<-G
# id<-rep(0,n)
# c<-1
# if(is.dag(g)){
#   while(!identical(n_p0,rep(-1,n))){
#     for (i in which(n_p0==0)) {
#       G0[i,]<-rep(0,p)
#       G0[i,i]<--1
#       id[i]<-c
#     }
#     c<-c+1
#     n_p0<-colSums(G0)
#   }
# }else{
#   print("This is not a DAG.")
# }
# topo<-order(id)
# 
# ## d.最大层级点
# max_id<-max(id)
# 
# ### topo_h:把前3层的节点序号换为0
# topo_h<-topo
# level_i_nodes <- sapply(1:3, function(i) which(id == i))
# topo_h[topo_h %in% unlist(level_i_nodes)] <- 0
# 
# 
# ### topo_no1: 节点序号按照拓扑排序排列的基础上，根节点为0
# topo_no1<-topo
# level_1_nodes <- which(id == 1)
# topo_no1[topo_no1 %in% level_1_nodes] <- 0
# 
# 
# ### topo_1：非根节点为0
# topo_1<-topo
# topo_1<-topo-topo_no1
