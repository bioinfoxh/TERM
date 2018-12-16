

# select candidates for seeds with various constrains: 
# cur_node: the current nodes of the n-module
# threshold: the values the candidate must satisfying 
# the strategy: 1. get the submatrix; 2. select only gene only one time based on the submatrix;


# networks=multi_netwk
# cur_node=ori_module[[i]]$members
# threshold=beta
# min_num=min_num_netw
# size_cand=size_cand


Cand_neibor <- function(networks,cur_node,threshold,min_num,size_cand)
{
  
  neibor=c();
  #print("yes, in the candidate");
  num_node = dim(networks)[1]
  num_netw = dim(networks)[3]
  rem_node  = setdiff(c(1:num_node),cur_node);
  num_rem_node = length(rem_node);
  rem_networks = networks[rem_node,cur_node,]; 
  #print(dim(rem_networks));
  rem_connect = array(0,dim=c(num_rem_node,num_netw+1));
  
  for(i in 1:num_netw){
    #print(  paste( "i for num_netw Cand_neibor: ",i,sep="")  )
    rem_connect[,i]=rowSums(rem_networks[,,i]);   
  }
  rem_connect[,num_netw+1]=rowSums(rem_connect[,1:num_netw])/num_netw;
  
  #******************* First level candidate ***********************************************************
  
  for(i in 1:num_rem_node){
    #print(  paste( "i for num_rem_node Cand_neibor: ",i,sep="")  )
    num_count = 0;
    index_count = rep(0,2);
    #for(tempi in 1:num_netw){ 
    if((max(rem_connect[i,1:(num_netw+1)])>threshold)&(min(rem_connect[i,1:(num_netw+1)])>0)){ 
      #num_count = num_count+1;
      num_count = num_count+sum(rem_connect[i,1:(num_netw+1)]>threshold)
    }  
    #}
    if(num_count>=min_num) {
      neibor=c(neibor,i); 
    }
  }
  
  #print(length(neibor));
  #******************* Second level candidate *********************************************************** 
  
  if( length(neibor)==1 ){
    neibor = rem_node[neibor]
  }
  
  if(length(neibor)>1 ){
    
    rem_node = rem_node[neibor];
    rem_networks = rem_networks[neibor,,];
    can_weight = rowSums(rem_networks[,,1]);
    for(tempi in 2:num_netw){ 
      can_weight = can_weight+rowSums(rem_networks[,,tempi]);  
    } 
    can_weight_list = order(can_weight,decreasing=TRUE);
    
    if(length(neibor)>size_cand) {
      neibor = rem_node[can_weight_list[1:size_cand]];
    }else{
      neibor = rem_node[can_weight_list];
    }
    
  }
  
  neibor;
}








# this is a fast version of N-modules: 
#cmatrix: the adjacent matrix, 
#cnode: the set of nodes; 
#curentropy: the current entropy for members 

# cmatrix=c_matrix
# cnode=c_node
# curentropy=ori_module[[i]]$p_entropy

multip_entropy <- function(cmatrix, cnode, curentropy, degree){
  
  #*************Step 1: compute the individual partial entropy for the candidate genes ********************
  #print("hello welcome to the fast one"); 
  num_mem = dim(curentropy)[1];                      # number of members
  num_can = length(cnode) - num_mem;                 # number of candidates
  num_netw = ncol(degree)
  
  #*******************Step 2: compute the  entropy changes for each input node ***************************************
  nsize = num_mem+1;
  Min_evalue = 1000000;
  Min_ematrix = c();
  Add_node = c();
  for(tempi in 1:num_can){        # for testing
    #print(  paste( "tempi for num_can multi_entropy: ",tempi,sep="")  )
    can_pentropy = array(0,dim=c(nsize,1+3*num_netw));
    
    can_pentropy[1:num_mem,] = curentropy;
    can_pentropy[nsize,1] = cnode[num_mem+tempi];
    can_pentropy[nsize,c(1:num_netw)+1+num_netw] = degree[can_pentropy[nsize,1],];
    #print(can_pentropy);
    rel_matrix = cmatrix[num_mem+tempi,,];
    can_pentropy[nsize,c(1:num_netw)+1] = colSums(rel_matrix);
    #print(rel_matrix);
    can_pentropy[1:num_mem,c(1:num_netw)+1] = can_pentropy[1:num_mem,c(1:num_netw)+1]+rel_matrix;  ### add the within degree of the new candidate to module members
    #print(can_pentropy);
    
    for(tempi1 in 1:nsize){
      for(tempi2 in 1:num_netw){
        
        temp_indegree = can_pentropy[tempi1,tempi2+1]
        temp_totaldegree = can_pentropy[tempi1,tempi2+1+num_netw]
        tempprob = temp_indegree/temp_totaldegree
        
        if(tempprob==0){
          can_pentropy[tempi1,tempi2+1+2*num_netw] = 2
        }else if( tempprob>=1 ){
          can_pentropy[tempi1,tempi2+1+2*num_netw] = 0
        }else if( (tempprob>0) & (tempprob<=0.5) ){
          can_pentropy[tempi1,tempi2+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob)
        }else if( (tempprob>0.5) & (tempprob<1)  ){
          can_pentropy[tempi1,tempi2+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob);  
        }
        
        
        #  if ( can_pentropy[tempi1,tempi2+1]==0){  
        #    can_pentropy[tempi1,tempi2+1+2*num_netw] = 2;  
        #    }
        # if( (can_pentropy[tempi1,tempi2+1]==can_pentropy[tempi1,tempi2+1+num_netw])&(can_pentropy[tempi1,tempi2+1]>0) ){
        #    can_pentropy[tempi1,tempi2+1+2*num_netw] = 0
        #    } 
        # if ( (2*can_pentropy[tempi1,tempi2+1]<=can_pentropy[tempi1,tempi2+1+num_netw])&(can_pentropy[tempi1,tempi2+1]>0) ){
        #    tempprob =  can_pentropy[tempi1,tempi2+1]/can_pentropy[tempi1,tempi2+1+num_netw];
        #    can_pentropy[tempi1,tempi2+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob);  
        #    } 
        # if( (2*can_pentropy[tempi1,tempi2+1]>can_pentropy[tempi1,tempi2+1+num_netw])&(can_pentropy[tempi1,tempi2+1]>0)&(can_pentropy[tempi1,tempi2+1]!=can_pentropy[tempi1,tempi2+1+num_netw]) ){  
        #    tempprob =  can_pentropy[tempi1,tempi2+1]/can_pentropy[tempi1,tempi2+1+num_netw]; 
        #    print(can_pentropy[tempi1,tempi2+1])
        #    print(can_pentropy[tempi1,tempi2+1+num_netw])
        #    print(paste( "tempprob greater than 0.5:", tempprob, sep="")  )
        #    can_pentropy[tempi1,tempi2+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob);  
        #    }
        
      }    # tempi2
    }   # tempi1 
    
    temp_entropy_value = sum(can_pentropy[,c(1:num_netw)+1+2*num_netw])/num_netw/(num_mem+1);
    temp_entropy_value[is.na(temp_entropy_value)] = 100;
    
    ori_entropy_value = colSums(curentropy[,c(1:num_netw)+1+2*num_netw])/num_mem;
    cur_entropy_value = colSums(can_pentropy[,c(1:num_netw)+1+2*num_netw])/(num_mem+1);
    threshold = min(ori_entropy_value-cur_entropy_value);
    if(length(cur_entropy_value[is.na(cur_entropy_value)])==0){
      if ((temp_entropy_value<Min_evalue)&(threshold>0.001)) { 
        Min_evalue = temp_entropy_value; 
        Min_ematrix = can_pentropy;  
        Add_node = cnode[tempi+num_mem]; 
      }
    }
  }
  #print(Add_node); 
  list(Min_evalue,Min_ematrix,Add_node);  
}




# the entropy value of each module

entropy<-function(module,networks){
  
  num_module = length(module);
  num_netw = dim(networks)[3];
  num_node = dim(networks)[1]
  
  for(i in 1:num_netw){ 
    diag(networks[,,i])=0 
  }
  
  degree=array(0,dim=c(num_node,num_netw));
  for(tempi in 1:num_netw){
    diag(networks[,,tempi])=0; 
    degree[,tempi]=rowSums(networks[,,tempi]); 
  }
  
  mentropy = array(0,dim=c(num_module,num_netw+1));
  wdegree = degree;
  for(i in 1:length(module)){
    #print(i);
    adjmatrix = networks[module[[i]]$members,module[[i]]$members,];
    modulelength = length(module[[i]]$members);
    for(j in 1:num_netw){
      tempmatrix = c();
      tempmatrix = adjmatrix[,,j];
      degree1 = wdegree[module[[i]]$members,j];#diag(tempmatrix); 
      tempmatrix[is.na(tempmatrix)] = 0;
      diag(tempmatrix) = 0;
      indegree = rowSums(tempmatrix);
      #print(degree);
      entropyvalue_k = c();
      for(k in 1:modulelength){
        
        #if(degree[k]==0){prob=0; print(tempmatrix);}else{prob = indegree[k]/degree[k];}
        prob = indegree[k]/degree1[k];
        if(prob==0){ 
          entropyvalue=2 
        }else if(prob>=1){
          entropyvalue=0
        }else if( (prob>0) & (prob<=0.5) ){   
          entropyvalue= 2+prob*log2(prob)+(1-prob)*log2(1-prob) 
        }else if( (prob>0.5) & (prob<1) ){
          entropyvalue = -prob*log2(prob)-(1-prob)*log2(1-prob)    
        }
        entropyvalue_k = c(entropyvalue_k,entropyvalue);
      }
      
      #mentropy[i,j]=log2(sum(entropyvalue_k)/modulelength);
      mentropy[i,j]=sum(entropyvalue_k)/modulelength;
      
    }
    mentropy[i,num_netw+1]=sum(mentropy[i,1:num_netw])/num_netw;
  }
  mentropy;
} 







# 7 = 1+3*num_netw
# 5 = 1+2*num_netw
# 3 = 1+num_netw
# c(2:3) = c(1:num_netw)+1
# c(4:5) = c(1:num_netw)+1+num_netw
# c(6:7) = c(1:num_netw)+1+2*num_netw

# the goal of the alogrithm is the n-modules in all these networkr 
nModule <- function(networks,seed_v){
  
  #options(warn=1)
  
  print("starting the n-module extraction precedure")
  adjmatrix = networks
  adjmatrix[is.na(adjmatrix)] = 0
  num_vertex = dim(networks)[1]
  num_netw = dim(networks)[3]
  num_node = dim(networks)[1];
  
  for(i in 1:num_netw){ 
    diag(networks[,,i])=0 
  }
  
  degree=array(0,dim=c(num_node,num_netw));
  for(tempi in 1:num_netw){
    diag(networks[,,tempi])=0; 
    degree[,tempi]=rowSums(networks[,,tempi]); 
  }
  
  ori_module <- vector(mode='list', length=length(seed_v)) # to store the n-modules
  print("finishing the module initail construction")
  
  #**************** initial procedure ***************************
  initial_index = 1  # 1: maximal 2: unoverlapping
  ori_module <- vector(mode='list', length=length(seed_v)) # to store the n-modules
  
  if(initial_index==1){
    
    for (i in 1:length(seed_v)){ 
      ori_module[[i]]<- list(name=paste("module", i, sep=" "), entropy=100,members=seed_v[i])
      tempmatrix = networks[ori_module[[i]]$members,,]  # co-expression weights of seed_v[i] in all the layers of networks
      
      # print(dim(tempmatrix)) 7737* 2
      can_weight = rowSums(tempmatrix)
      if( sum(can_weight)==0 ){ next }
      
      can_gene_list = order(can_weight,decreasing=TRUE)
      ori_module[[i]]$members=c(ori_module[[i]]$members, can_gene_list[1])
      
      tempmatrix = networks[,ori_module[[i]]$members,]
      can_weight = rowSums(tempmatrix[,,1])
      
      for(tempi in 2:num_netw){ 
        can_weight = can_weight+rowSums(tempmatrix[,,tempi])  
      }
      if( sum(can_weight)==0 ){ next }
      
      can_gene_list = order(can_weight,decreasing=TRUE)
      
      for(tempi in 1:length(can_gene_list)){
        if (!(can_gene_list[tempi] %in% ori_module[[i]]$members)){  
          ori_module[[i]]$members=c(ori_module[[i]]$members, can_gene_list[tempi])  
          break    ### add the one with the highest weight into module i, and then get out of the loop
        }  
      } 
      
      #  print(ori_module[[i]]$members)
      
      ori_module[[i]]$members = sort( ori_module[[i]]$members)
      #m_entropy=array(0,dim=c(length( ori_module[[i]]$members),7)) #c1 genes,c2-c3 within-module degree, c4-c5 total degree, c6-c7 entropy
      ##### !!!!!!!! WARNING: 7 is for two-layer network, for more than 2 layers, it should be modified
      m_entropy=array(0,dim=c(length( ori_module[[i]]$members), 1+3*num_netw)) #### for multi-layer networks with more than 2 layers
      m_entropy[,1]=ori_module[[i]]$members                        #c1: genes
      m_entropy[,c(1:num_netw)+1+num_netw]=degree[ori_module[[i]]$members,]          #total degrees in each layer 
      
      for(tempj in 1:num_netw){
        m_entropy[,tempj+1]=rowSums(networks[ori_module[[i]]$members,ori_module[[i]]$members,tempj]) #within module degrees in each layer
        
        # compute the entropy for each gene
        for(tempj1 in 1:dim(m_entropy)[1]){
          
          temp_indegree = m_entropy[tempj1,tempj+1]
          temp_totaldegree = m_entropy[tempj1,tempj+1+num_netw]
          tempprob = temp_indegree/temp_totaldegree
          
          if(tempprob==0){
            m_entropy[tempj1,tempj+1+2*num_netw] = 2
          }else if( tempprob>=1 ){
            m_entropy[tempj1,tempj+1+2*num_netw] = 0
          }else if( (tempprob>0) & (tempprob<=0.5) ){
            m_entropy[tempj1,tempj+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob)
          }else if( (tempprob>0.5) & (tempprob<1)  ){
            m_entropy[tempj1,tempj+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob)  
          }
          
          
          # if ( m_entropy[tempj1,tempj+1] == 0 ){ 
          #   m_entropy[tempj1,tempj+1+2*num_netw] = 2  
          #   }else if ( (m_entropy[tempj1,tempj+1]==m_entropy[tempj1,tempj+1+num_netw]) ){
          #       m_entropy[tempj1,tempj+1+2*num_netw] = 0
          #     }else if ( (2*m_entropy[tempj1,tempj+1] <= m_entropy[tempj1,tempj+1+num_netw]) ) { 
          #          tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+1+num_netw] 
          #          m_entropy[tempj1,tempj+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob)  
          #       }else{  
          #         tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+1+num_netw] 
          #         m_entropy[tempj1,tempj+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob)   
          #         }
          
          
        }#tempj1
      } #tempj
      #print(ori_module[[i]]$members)
      ori_module[[i]]$p_entropy = m_entropy
      # print(ori_module[[i]]$p_entropy)
      ori_module[[i]]$v_entropy = sum(m_entropy[,c(1:num_netw)+1+2*num_netw])/num_netw/length(ori_module[[i]]$members)
    }
    
  }else{
    
    init_index_vect = rep(0,num_vertex)
    
    for (i in 1:length(seed_v)){ 
      ori_module[[i]]<- list(name=paste("module", i, sep=" "), entropy=100,members=seed_v[i])
      tempmatrix = networks[ori_module[[i]]$members,,]
      
      # print(dim(tempmatrix)) 7737* 2
      can_weight = rowSums(tempmatrix)
      if( sum(can_weight)==0 ){ next }
      
      can_gene_list = order(can_weight,decreasing=TRUE)
      
      temp.index = 1
      while(length(ori_module[[i]]$members)<3){  
        if(!init_index_vect[can_gene_list[temp.index]] ){  
          if( can_weight[ can_gene_list[temp.index] ]==0   ){
            ori_module[[i]]$members = c(ori_module[[i]]$members,NULL)
            break
          }else{
            ori_module[[i]]$members = c(ori_module[[i]]$members,can_gene_list[temp.index]) 
            init_index_vect[can_gene_list[temp.index]]=1 
          }
        }
        temp.index = temp.index +1
      }     
      
      if(length(ori_module[[i]]$members)<3){next}
      
      ori_module[[i]]$members = sort( ori_module[[i]]$members)
      m_entropy=array(0,dim=c(length( ori_module[[i]]$members),1+3*num_netw)) #c1:genesc2-c3: in degree, c4-c5:  degree, c6-c7:entropy
      m_entropy[,1]=ori_module[[i]]$members                        #c1: genes
      m_entropy[,c(1:num_netw)+1+num_netw]=degree[ori_module[[i]]$members,]          #c2-c3: in degree
      for(tempj in 1:num_netw){
        m_entropy[,tempj+1]=rowSums(networks[ori_module[[i]]$members,ori_module[[i]]$members,tempj]) #c
        
        # compute the entropy for each gene
        for(tempj1 in 1:dim(m_entropy)[1]){
          
          temp_indegree = m_entropy[tempj1,tempj+1]
          temp_totaldegree = m_entropy[tempj1,tempj+1+num_netw]
          tempprob = temp_indegree/temp_totaldegree
          
          if(tempprob==0){
            m_entropy[tempj1,tempj+1+2*num_netw] = 2
          }else if( tempprob>=1 ){
            m_entropy[tempj1,tempj+1+2*num_netw] = 0
          }else if( (tempprob>0) & (tempprob<=0.5) ){
            m_entropy[tempj1,tempj+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob)
          }else if( (tempprob>0.5) & (tempprob<1)  ){
            m_entropy[tempj1,tempj+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob) 
          }
          
          # if (m_entropy[tempj1,tempj+1] == 0){ 
          #   m_entropy[tempj1,tempj+1+2*num_netw] = 2  
          #   }else if( (m_entropy[tempj1,tempj+1]==m_entropy[tempj1,tempj+1+num_netw]) ){
          #        m_entropy[tempj1,tempj+1+2*num_netw] = 0
          #      }else if ( (2*m_entropy[tempj1,tempj+1] <= m_entropy[tempj1,tempj+1+num_netw]) ){ 
          #          tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+1+num_netw] 
          #          m_entropy[tempj1,tempj+1+2*num_netw] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob)  
          #        }else{  
          #          tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+1+num_netw]
          #          m_entropy[tempj1,tempj+1+2*num_netw] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob)   
          #       }
          #     
          
        } #tempj1
      } #tempj
      #print(ori_module[[i]]$members)
      ori_module[[i]]$p_entropy = m_entropy
      # print(ori_module[[i]]$p_entropy)
      ori_module[[i]]$v_entropy = sum(m_entropy[,c(1:num_netw)+1+2*num_netw])/num_netw/length(ori_module[[i]]$members)
    }
    
  }
  #*************************************************************
  #***** The candidates should meet two requirements:
  #*****   1. should be highly co-expressed 
  #*****   2. sensitivity to the drug 
  
  beta = 0.01
  alpha = 0.1         # the threshold for the 
  min_num_netw = 1
  size_cand = 100
  beta2 = 0.001
  m_size = 200 # the size of the 
  
  #**************** expand procedure ***************************
  
  for (i in 1:length(seed_v)){ #for checking
    
    # tempindex indicates the maximum size of a module 
    
    print(paste("extracting the module ", i," / ",length(seed_v), sep=""))
    
    if(length(ori_module[[i]]$members)<3){next}
    
    tempindex=0
    
    #temp_loop_count=0
    while((length(ori_module[[i]]$members)<m_size)&(tempindex<1)){
      #temp_loop_count = temp_loop_count+1  ##### for test
      #print(paste("temp_loop_count= ",temp_loop_count,sep="" )  )
      
      #source("./Cand_neibor.R")
      neibor = Cand_neibor(networks,ori_module[[i]]$members,beta,min_num_netw,size_cand) 
      
      if (length(neibor)==0){ 
        print("no candidates") 
        tempindex=101  
      }else{
        c_node = c(ori_module[[i]]$members,neibor)
        c_matrix = networks[c_node,ori_module[[i]]$members,]
        #source("./multip_entropy.R")
        candid =multip_entropy(c_matrix,c_node,ori_module[[i]]$p_entropy,degree) #c_node is the union of $members and neibors
        
        #print(candid[[1]])
        if ((ori_module[[i]]$v_entropy-candid[[1]])<beta2|is.na(ori_module[[i]]$v_entropy-candid[[1]])){ 
          #print("the entropy is: ")
          #print(ori_module[[i]]$v_entropy)
          #print(candid[[1]])
          #print("the new node")
          # print(candid[[3]]) 
          #  print("no further improvement")  
          tempindex=101 
        }else{  
          # print("the entropy is: ")
          #print(ori_module[[i]]$v_entropy)
          #print(candid[[1]])
          # print("the new node")
          #  print(candid[[3]])
          ori_module[[i]]$members = c(ori_module[[i]]$members,candid[[3]])
          ori_module[[i]]$members = sort( ori_module[[i]]$members)
          tempmatrix = c()
          tempmatrix = candid[[2]]
          ori_module[[i]]$members=sort(tempmatrix[,1])
          tempmatrix= tempmatrix[order(tempmatrix[,1]),]
          ori_module[[i]]$p_entropy = c()
          ori_module[[i]]$p_entropy = tempmatrix
          ori_module[[i]]$v_entropy=candid[[1]]
        }  
        
      }
      
    }
    
    ori_module[[i]]$matrix = networks[sort(ori_module[[i]]$members),sort(ori_module[[i]]$members),]
    ori_module[[i]]$entropy = colSums(ori_module[[i]]$p_entropy[,c(1:num_netw)+1+2*num_netw])/length(ori_module[[i]]$members)
    
   # print(ori_module[[i]]$members)
    
  }
  #************************************************************* 
  ori_module
}






