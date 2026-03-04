
rm(list = ls())


### 1.??????
library(rstan)
library(ggplot2)
library(loo)
options(mc.cores = (parallel::detectCores())-2)
rstan_options(auto_write = TRUE)
library(igraph)
library(MASS)
library(openxlsx)

### 2.??????s
source('gsp_function/random_graph.R')
source('gsp_function/convert_to_moral_graph.R')
source('gsp_function/topo_order.R')


filename<-'new_simulation1_graph_rep20_noise1'
dir.create(filename)
dir.create(paste(filename,'/data',sep=''))
dir.create(paste(filename,'/rstan_result',sep=''))
dir.create(paste(filename,'/data_gsp',sep=''))
dir.create(paste(filename,'/mse_result',sep=''))

# n_list<-c(10,30,50)
n_list<-c(50)
# m_list<-rep(n_list,each=3)*rep(c(1,2,4),3)
m_list <- c(150)
J<-1
K<-1
seed<-129
nsim<-3*3*20

para_list<-c('tau','psi',"eta","w","sigma2","phi","h","Gdiff_1","Gdiff_2","mu")
#para_list<-c('tau','psi',"delta","eta","w","sigma2","phi","h","Gdiff_1","Gdiff_2","mu")
#para_list<-c("delta","w","sigma2","h","Gdiff_1","Gdiff_2","mu")

rep_times<-2
colname<-c('n','m','ini_sigma','node_sigma','noise_sigma','rep'
           ,'X_mse','gsp_mse','x_snr','gsp_snr')
#???ڵ???׼??
ini_sigma<-10
#?ڵ??Ŷ?
node_sigma<-1
#????
noise_sigma<-1
mse_result<-cbind(n=rep(n_list,each=3*rep_times)
                  ,m=rep(m_list,each=rep_times)
                  ,ini_sigma=rep(ini_sigma,nsim)
                  ,node_sigma=rep(node_sigma,nsim)
                  ,noise_sigma=rep(noise_sigma,nsim)
                  ,rep=rep(c(1:rep_times),nsim/rep_times))

#????????ͼ
gsp_mse<-NULL
x_mse<-NULL
gsp_snr<-NULL
x_snr<-NULL


for(n in n_list){
  m_list<-c(n,2*n,4*n)
  for(m in m_list){
    for(rep in 1:rep_times){
      time0 <- Sys.time()
      seed<-10*(n+m+rep)
      set.seed(seed)
      G<-random_graph(n=n,m=m)
      g<-graph_from_adjacency_matrix(G,mode = 'directed')
      #?????ڽӾ???
      write.csv(G,paste0(filename,"/G_sim1_n",n,"_m",m,"_rep",rep,".csv"),row.names = F)
      #?????��?ͼ
      MoG<-unname(as.matrix(G))
      MoG<-convert_to_moral_graph(MoG)
      write.csv(MoG,paste0(filename,"/MoG_sim1_n",n,"_m",m,"_rep",rep,".csv"),row.names = F)
      ############
      source('gsp_function/random_graph.R')
      source('gsp_function/convert_to_moral_graph.R')
      source('gsp_function/topo_order.R')
      ############
      ## ????
      topo_order<-topo_order(g,G,n,p=1,h=2)
      id<-topo_order$id
      topo<-topo_order$topo
      n_p<-topo_order$n_p
      topo_1<-topo_order$topo_1
      topo_no1<-topo_order$topo_no1
      topo_h<-topo_order$topo_h
      topo_2<-topo_order$topo_2
      ############
      set.seed(seed)
      inital<-rnorm(length(which(id==1)),0,ini_sigma)
      names(inital)<-which(id==1)
      mu<-matrix(0,J,n)
      for (i in topo) {
        if(i %in% which(id==1)){
          set.seed(seed+i)
          mu[,i]<-unname(inital[as.character(i)])+rnorm(J,0,node_sigma)
        }else{
          set.seed(seed+i)
          mu[,i]<-mu[,]%*%G[,i]/n_p[i]+rnorm(J,0,node_sigma)
          
        }
      }
      set.seed(seed)
      X<-matrix(rnorm(J*n,0,noise_sigma),J,n)
      for (i in 1:n) {
        X[,i]<-mu[,i]+X[,i]
      }
      ## save data
      data<-data.frame()
      X0<-t(X)
      mu0<-t(mu)
      data<-data.frame(node=rep(1:n),X=as.vector(X0),mu=as.vector(mu0))
      write.csv(data,paste0(filename,paste('/data/data_n',n,"_m",m,"_rep",rep,sep=''),'.csv'),row.names = F)
      cauchy_par<-sd(X)/sqrt(n)
      data_GSP<-list(n=n,J=J,K=K,X=X,G=G,n_p=n_p,topo=topo,topo_h=topo_h,topo_2=topo_2,topo_1=topo_1,topo_no1=topo_no1,cauchy_par=cauchy_par)
      GSP<-stan_model("Z:/U_nodn")
      
      start_time <- Sys.time() 
      fit<-sampling(GSP,data_GSP,chain=4,iter=2000)
      end_time <- Sys.time() 
      print(end_time - start_time)
      ### 5.????????????
      wb <- createWorkbook()
      for(w in 1:length(para_list)){
        addWorksheet(wb, sheetName =para_list[w])
        writeData(wb, sheet =w, summary(fit,par=para_list[w])$summary)
      }
      saveWorkbook(wb, paste0(filename,paste('/rstan_result/rstan_result_n',n,"_m",m,"_rep",rep,sep=''),'.xlsx'), overwrite = TRUE)
      
      ### 6.????mu
      rstan_mu<-summary(fit,par='mu')$summary[,1]
      rstan_mu<-matrix(rstan_mu,ncol=1)
      rstan_mu<-as.vector(rstan_mu)
      data_gsp<-cbind(data,mu_GSP=rstan_mu)
      write.csv(data_gsp,paste0(filename,paste('/data_gsp/data_gsp_n',n,"_m",m,"_rep",rep,sep=''),'.csv'),row.names = F)
      ##????mse????
      rep_gsp_mse<-mean((data_gsp$mu-data_gsp$mu_GSP)^2)
      rep_x_mse<-mean((data_gsp$mu-data_gsp$X)^2)
      rep_gsp_snr<-10*log10(sum((data_gsp$mu)^2)/sum((data_gsp$mu-data_gsp$mu_GSP)^2))
      rep_x_snr<-10*log10(sum((data_gsp$mu)^2)/sum((data_gsp$mu-data_gsp$X)^2))
      gsp_mse<-c(gsp_mse,rep_gsp_mse)
      x_mse<-c(x_mse,rep_x_mse)
      print(x_mse)
      gsp_snr<-c(gsp_snr,rep_gsp_snr)
      x_snr<-c(x_snr,rep_x_snr)
    }
    
      run_time <- Sys.time() - time0
      
  }
}
mse_result<-cbind(mse_result,x_mse,x_snr,gsp_mse,gsp_snr)
colnames(mse_result)<-colname
write.csv(mse_result,paste0(filename,paste('/mse_result/mse_result.csv')),row.names = F)
