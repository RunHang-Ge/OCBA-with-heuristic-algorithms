################# 新OCBA连续蚁群算法测试 ##################
elite_new_partrank_ttest_OCBA = function(a_test){
  #calcu_iter = para
  #max_iter = 100
  #temp_six_avg = matrix(data=NA,nrow = max_iter,ncol = calcu_iter)
  #temp_six_true = matrix(data=NA,nrow = max_iter,ncol = calcu_iter)
  #for (a in 1:calcu_iter){
  ### 定义参数
  k=50 #k为群体(design)数量
  dimension=10 #n为解的维数
  zeta = 0.85 #ζ为方差收敛参数
  q = 0.1 #q为权重参数
  m = 10 #m为每次更新的数量
  lim = 13 #lim为全排列的阈值
  max_iter = 4000 #最大迭代次数
  cf = 2.5 #精英选择的常数
  ci = 0.5 #精英选择的常数
  n0 = 10 #初始每个design的计算量
  n1 = 100 #最后筛选时的计算量
  n = rep(n0,k+1) #储存各个design的计算量
  history_n = matrix(data=NA,ncol = max_iter, nrow = k+1)
  N = 5000 #计算量分配阈值
  delta = 100 #每次迭代增加计算量Δ
  sigma = 10 #计算误差大小
  object_func = function(x_one,x_two,x_three,x_four,x_five,x_six,x_seven,x_eight,x_nine,x_ten,sigma_test){ #目标函数
    #temp_one = (1+(x_one+x_two+1)^2*(19-14*x_one+3*x_one^2-14*x_two+6*x_one*x_two+3*x_two^2))*(30+(2*x_one-3*x_two)^2*(18-32*x_one+12*x_one^2+48*x_two-36*x_one*x_two+27*x_two^2))+rnorm(1,mean = 0,sd = sigma_test)  #Goldstein and Price 效果不错
    #temp_one = x_one^2+2*x_two^2-(3/10)*cos(3*pi*x_one)-(2/5)*cos(4*pi*x_two)+7/10+rnorm(1,mean = 0,sd = sigma_test) #B2 效果不好
    #temp_one = (x_two-(5/(4*pi^2)*x_one^2)+5/pi*x_one-6)^2+10*(1-(1/(8*pi)))*cos(x_one)+10+rnorm(1,mean = 0,sd = sigma_test) #Branin 效果不好，趋于最优值速度慢
    #temp_one = 4*x_one-2.1*x_one^4+1/3*x_one^6+x_one*x_two-4*x_two^2+4*x_two^4+1.0316+rnorm(1,mean = 0,sd = sigma_test)  #Six-hump 效果不错
    #temp_one = -cos(x_one)*cos(x_two)*exp(-((x_one-pi)^2+(x_two-pi)^2))+rnorm(1,mean = 0,sd = sigma_test) #Easom 效果不好
    #temp_one = (x_one-x_two)^2+((x_one+x_two-10)/3)^2+rnorm(1,mean = 0,sd = sigma_test) #Martin and gaddy
    #temp_one = x_one^2+x_two^2+rnorm(1,mean = 0,sd = sigma_test)
    
    temp_one = -20*exp(-0.2*sqrt((x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2)/10))-exp((cos(2*pi*x_one)+cos(2*pi*x_two)+cos(2*pi*x_three)+cos(2*pi*x_four)+cos(2*pi*x_five)+cos(2*pi*x_six)+cos(2*pi*x_seven)+cos(2*pi*x_eight)+cos(2*pi*x_nine)+cos(2*pi*x_ten))/10)+20+exp(1)+rnorm(1,mean = 0,sd = sigma_test) #Ackley
    #temp_one = 100*((x_one^2-x_two)^2+(x_two^2-x_three)^2+(x_three^2-x_four)^2+(x_four^2-x_five)^2+(x_five^2-x_six)^2+(x_six^2-x_seven)^2+(x_seven^2-x_eight)^2+(x_eight^2-x_nine)^2+(x_nine^2-x_ten)^2)+(x_one-1)^2+(x_two-1)^2+(x_three-1)^2+(x_four-1)^2+(x_five-1)^2+(x_six-1)^2+(x_seven-1)^2+(x_eight-1)^2+(x_nine-1)^2+(x_ten-1)^2+rnorm(1,mean = 0,sd = sigma_test) #Rosenbrock
    #temp_one = sum(abs(x_one)+abs(x_two)+abs(x_three)+abs(x_four)+abs(x_five)+abs(x_six)+abs(x_seven)+abs(x_eight)+abs(x_nine)+abs(x_ten))+abs(x_one)*abs(x_two)*abs(x_three)*abs(x_four)*abs(x_five)*abs(x_six)*abs(x_seven)*abs(x_eight)*abs(x_nine)*abs(x_ten)+rnorm(1,mean = 0,sd = sigma_test)
    #temp_one = (x_one^2-10*cos(2*pi*x_one)+10)+(x_two^2-10*cos(2*pi*x_two)+10)+(x_three^2-10*cos(2*pi*x_three)+10)+(x_four^2-10*cos(2*pi*x_four)+10)+(x_five^2-10*cos(2*pi*x_five)+10)+(x_six^2-10*cos(2*pi*x_six)+10)+(x_seven^2-10*cos(2*pi*x_seven)+10)+(x_eight^2-10*cos(2*pi*x_eight)+10)+(x_nine^2-10*cos(2*pi*x_nine)+10)+(x_ten^2-10*cos(2*pi*x_ten)+10)+rnorm(1,mean = 0,sd = sigma_test) #Rastrigin
    #temp_one = x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2+rnorm(1,mean = 0,sd = sigma_test) #Sphere
    #temp_one = x_one^2+(1e02^(1/9)*x_two)^2+(1e02^(2/9)*x_three)^2+(1e02^(3/9)*x_four)^2+(1e02^(4/9)*x_five)^2+(1e02^(5/9)*x_six)^2+(1e02^(6/9)*x_seven)^2+(1e02^(7/9)*x_eight)^2+(1e02^(8/9)*x_nine)^2+(1e02*x_ten)^2+rnorm(1,mean = 0,sd = sigma_test)# shifted Rotated
    #temp_one = 1/4000*(x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2)-(cos(x_one/sqrt(1))*cos(x_two/sqrt(2))*cos(x_three/sqrt(3))*cos(x_four/sqrt(4))*cos(x_five/sqrt(5))*cos(x_six/sqrt(6))*cos(x_seven/sqrt(7))*cos(x_eight/sqrt(8))*cos(x_nine/sqrt(9))*cos(x_ten/sqrt(10)))+1+rnorm(1,mean = 0,sd = sigma_test)
    return(temp_one)
  }
  object_func_true = function(x_one,x_two,x_three,x_four,x_five,x_six,x_seven,x_eight,x_nine,x_ten,sigma_test){ #目标函数真值
    #temp_six = (x_two-(5/(4*pi^2)*x_one^2)+5/pi*x_one-6)^2+10*(1-(1/(8*pi)))*cos(x_one)+10 #Branin
    #temp_six = (1+(x_one+x_two+1)^2*(19-14*x_one+3*x_one^2-14*x_two+6*x_one*x_two+3*x_two^2))*(30+(2*x_one-3*x_two)^2*(18-32*x_one+12*x_one^2+48*x_two-36*x_one*x_two+27*x_two^2))
    #temp_six = 4*x_one-2.1*x_one^4+1/3*x_one^6+x_one*x_two-4*x_two^2+4*x_two^4+1.0316  #Six-hump
    #temp_six = x_one^2+2*x_two^2-(3/10)*cos(3*pi*x_one)-(2/5)*cos(4*pi*x_two)+7/10 #B2
    #temp_six = -cos(x_one)*cos(x_two)*exp(-((x_one-pi)^2+(x_two-pi)^2))
    #temp_six = (x_one-x_two)^2+((x_one+x_two-10)/3)^2
    #temp_six = x_one^2+x_two^2
    
    #temp_six = x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2
    temp_six = -20*exp(-0.2*sqrt((x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2)/10))-exp((cos(2*pi*x_one)+cos(2*pi*x_two)+cos(2*pi*x_three)+cos(2*pi*x_four)+cos(2*pi*x_five)+cos(2*pi*x_six)+cos(2*pi*x_seven)+cos(2*pi*x_eight)+cos(2*pi*x_nine)+cos(2*pi*x_ten))/10)+20+exp(1) #Ackley
    #temp_six = 100*((x_one^2-x_two)^2+(x_two^2-x_three)^2+(x_three^2-x_four)^2+(x_four^2-x_five)^2+(x_five^2-x_six)^2+(x_six^2-x_seven)^2+(x_seven^2-x_eight)^2+(x_eight^2-x_nine)^2+(x_nine^2-x_ten)^2)+(x_one-1)^2+(x_two-1)^2+(x_three-1)^2+(x_four-1)^2+(x_five-1)^2+(x_six-1)^2+(x_seven-1)^2+(x_eight-1)^2+(x_nine-1)^2+(x_ten-1)^2 #Rosenbrock
    #temp_six = sum(abs(x_one)+abs(x_two)+abs(x_three)+abs(x_four)+abs(x_five)+abs(x_six)+abs(x_seven)+abs(x_eight)+abs(x_nine)+abs(x_ten))+abs(x_one)*abs(x_two)*abs(x_three)*abs(x_four)*abs(x_five)*abs(x_six)*abs(x_seven)*abs(x_eight)*abs(x_nine)*abs(x_ten)
    #temp_six = (x_one^2-10*cos(2*pi*x_one)+10)+(x_two^2-10*cos(2*pi*x_two)+10)+(x_three^2-10*cos(2*pi*x_three)+10)+(x_four^2-10*cos(2*pi*x_four)+10)+(x_five^2-10*cos(2*pi*x_five)+10)+(x_six^2-10*cos(2*pi*x_six)+10)+(x_seven^2-10*cos(2*pi*x_seven)+10)+(x_eight^2-10*cos(2*pi*x_eight)+10)+(x_nine^2-10*cos(2*pi*x_nine)+10)+(x_ten^2-10*cos(2*pi*x_ten)+10) #Rastrigin
    #temp_six = x_one^2+(1e02^(1/9)*x_two)^2+(1e02^(2/9)*x_three)^2+(1e02^(3/9)*x_four)^2+(1e02^(4/9)*x_five)^2+(1e02^(5/9)*x_six)^2+(1e02^(6/9)*x_seven)^2+(1e02^(7/9)*x_eight)^2+(1e02^(8/9)*x_nine)^2+(1e02*x_ten)^2
    #temp_six = 1/4000*(x_one^2+x_two^2+x_three^2+x_four^2+x_five^2+x_six^2+x_seven^2+x_eight^2+x_nine^2+x_ten^2)-(cos(x_one/sqrt(1))*cos(x_two/sqrt(2))*cos(x_three/sqrt(3))*cos(x_four/sqrt(4))*cos(x_five/sqrt(5))*cos(x_six/sqrt(6))*cos(x_seven/sqrt(7))*cos(x_eight/sqrt(8))*cos(x_nine/sqrt(9))*cos(x_ten/sqrt(10)))+1
    return(temp_six)
  }
  x_max = 10 #搜索上线
  x_min = -10 #搜索下限
  
  ### 均匀分布生成初始个体
  pop = matrix(data = NA,nrow = k,ncol = dimension)
  for (i in 1:dimension){
    pop[,i] = runif(k,min = x_min,max = x_max)
  }
  
  ### 计算每个design的μ和σ并排序
  pop_value = matrix(data = NA,nrow = k,ncol = 4) #4列分别储存序号、μ、σ、真值
  pop_value[,1] = seq(1,k)
  temp_two = matrix(data = NA,ncol = n0,nrow = k)
  for (i in 1:k){ #进行随机仿真
    for (j in 1:n0){
      temp_two[i,j] = object_func(pop[i,1],pop[i,2],pop[i,3],pop[i,4],pop[i,5],pop[i,6],pop[i,7],pop[i,8],pop[i,9],pop[i,10],sigma)
    }
  }
  for (i in 1:k){ #计算μ和σ,真值
    pop_value[i,2] = abs(mean(temp_two[i,]))
    pop_value[i,3] = sd(temp_two[i,])
    pop_value[i,4] = object_func_true(pop[i,1],pop[i,2],pop[i,3],pop[i,4],pop[i,5],pop[i,6],pop[i,7],pop[i,8],pop[i,9],pop[i,10],sigma)
  }
  pop_value_sorted = pop_value[order(pop_value[,2]),]
  
  ### 记录初代最优个体，μ，真值
  history_best = matrix(data = NA,nrow = max_iter,ncol = dimension+2)
  history_best[1,] = c(pop[pop_value_sorted[1,1],],pop_value_sorted[1,2],pop_value_sorted[1,4])
  
  
  ############# ACO迭代开始
  #ACO_OCBA = function(z_test){
  z = 2
  while (z <= max_iter){
    ############# OCBA开始
    ### 记录sorted后的计算量比例α,最后一个永远是上代最优的α
    alpha = rep(1/(k+1),k+1)
    
    ### 计算种群在计算量均为n0的排名、μ、σ，得到当代最优b。此处pop_value的最后一个值一定是上一代最优
    pop_value = matrix(data = NA,nrow = k+1,ncol = 4) #4列分别储存序号、μ、σ、真值
    pop_value[,1] = seq(1,k+1)
    temp_two = matrix(data = NA,ncol = n0,nrow = k+1) #所有仿真值的记录
    for (i in 1:(k+1)){ #进行随机仿真
      if (i == k+1){ #这个是上一代最优，第一次迭代的上次最优是history_best
        for (j in 1:n0){
          temp_two[i,j] = object_func(history_best[z-1,1],history_best[z-1,2],history_best[z-1,3],history_best[z-1,4],history_best[z-1,5],history_best[z-1,6],history_best[z-1,7],history_best[z-1,8],history_best[z-1,9],history_best[z-1,10],sigma)
        }
        pop_value[i,2] = abs(mean(temp_two[i,]))
        pop_value[i,3] = sd(temp_two[i,])
        pop_value[i,4] = object_func_true(history_best[z-1,1],history_best[z-1,2],history_best[z-1,3],history_best[z-1,4],history_best[z-1,5],history_best[z-1,6],history_best[z-1,7],history_best[z-1,8],history_best[z-1,9],history_best[z-1,10],sigma)
      }else{ #这些是目前的种群
        for (j in 1:n0){
          temp_two[i,j] = object_func(pop[i,1],pop[i,2],pop[i,3],pop[i,4],pop[i,5],pop[i,6],pop[i,7],pop[i,8],pop[i,9],pop[i,10],sigma)
        }
        pop_value[i,2] = abs(mean(temp_two[i,]))
        pop_value[i,3] = sd(temp_two[i,])
        pop_value[i,4] = object_func_true(pop[i,1],pop[i,2],pop[i,3],pop[i,4],pop[i,5],pop[i,6],pop[i,7],pop[i,8],pop[i,9],pop[i,10],sigma)
      }
    }
    pop_value_sorted = rbind(pop_value[order(pop_value[-(k+1),2]),],pop_value[k+1,]) #前k行以仿真值排序,最后一行仍是上代最优
    b = pop_value_sorted[1,1] #记录本代最优
    best_before = history_best[z-1,1:dimension] #记录上代最优
    
    ### 增加n的对应序号
    n = cbind(pop_value_sorted[,1],n)
    sum_N = sum(n[,2])
    
    ############# 初始t检验
    t_result = rep(NA,lim-1)
    df = rep(NA,lim-1)
    Q_set = matrix(data = NA, ncol=2,nrow = lim-1)
    for (i in 1:(lim-1)){
      t_result[i] = (pop_value_sorted[i+1,2]-pop_value_sorted[1,2])/sqrt(pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))
      df[i] = (pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/((pop_value_sorted[1,3]^2/(alpha[1]*sum_N))^2/(alpha[1]*sum_N-1)+(pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/(alpha[i+1]*sum_N-1))
      if (t_result[i] < qt(0.95,df = df[i])){
        Q_set[i,1] = pop_value_sorted[i+1,2]
        Q_set[i,2] = "Same"
      }else{
        Q_set[i,2] = "Not Same"
      }
    }
    Q_mean = sum(pop_value_sorted[c(1,which(Q_set[,2]=="Same")+1),2])/length(c(1,which(Q_set[,2]=="Same")+1)) ###更换均值
    ### 注意，这里的第一个是排名第二的个体与最优个体的对比
    
    
    ############# 迭代开始
    n_temp = rep(n0,k+1)
    while(sum_N < N){
      ### 计算阈值M
      M=(pop_value_sorted[lim+1,3]*pop_value_sorted[lim,2]/sqrt(alpha[lim+1]*sum_N)+pop_value_sorted[lim,3]*pop_value_sorted[lim+1,2]/sqrt(alpha[lim]*sum_N))/(pop_value_sorted[lim+1,3]/sqrt(alpha[lim+1]*sum_N)+pop_value_sorted[lim,3]/sqrt(alpha[lim]*sum_N))
      
      ### 判断SG是否为Φ
      if (pop_value_sorted[1,2] >= pop_value_sorted[k+1,2]){ #如果这一代最优不如上一代最优
        converge_rate = matrix(data = NA,nrow = k,ncol = 4) #用来保存三个值和类别
        for (i in 1:k){ #计算Ii+1,Ii-1,G
          if (i == 1){
            converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]+1e12)^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
          }else if (i <= lim-1){
            if (Q_set[i-1,2] == "Same"){
              converge_rate[i,1] = alpha[i]*(Q_mean-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(Q_mean-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (Q_mean-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            }else{
              converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,2]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            }
          }else if (i == lim){
            if (Q_set[i-1,2] == "Same"){
              converge_rate[i,1] = alpha[i]*(Q_mean-M)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(Q_mean-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (Q_mean-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            }else{
              converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-M)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            }
          }else{
            converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-M)^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,2] = (pop_value_sorted[i,2]-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            converge_rate[i,3] = 1e12
          }
        }
        
        sum_alpha = 0 #用于计算alpha[k+1]
        for (i in 1:k){ #先计算种群的alpha
          if (i == 1){
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
              converge_rate[i,4] = 1
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i+1,2]))^2
            }
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
              converge_rate[i,4] = 2
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i-1,2]))^2
            }
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
              converge_rate[i,4] = 3
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[k+1,2]))^2
              sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
            }
          }else if (i <= lim-1){
            if (Q_set[i-1,2] == "Same"){
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i+1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[k+1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }else{
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i+1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[k+1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }
          }else if(i == lim){
            if (Q_set[i-1,2] == "Same"){
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-M))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[k+1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }else{
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-M))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[k+1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }
          }else{
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
              converge_rate[i,4] = 1
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-M))^2
            }
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
              converge_rate[i,4] = 2
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[k+1,2]))^2
              sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
            }
          }
        }
        alpha[k+1] = pop_value_sorted[k+1,3]*sqrt(sum_alpha)
        
        ### 开始增加计算量,注意此时的α是sorted之后的
        n_new = sum(c(n[,2],delta))/sum(alpha)*alpha
        n_add = rep(NA,k+1)
        for (i in 1:(k+1)){
          if (n_new[i] <= n[i,2]){
            n_add[i] = 0
          }else{
            n_add[i] = n_new[i]-n[i,2]
          }
        }
        n_add = round(n_add/sum(n_add)*delta)
        #View(n_add)
        ### 新增计算量后的排名、μ、σ，得到当代最优b。此处pop_value的最后一个值一定是上一代最优
        ### 由于我们不需要更换上一代最优，因此直接增加仿真数量即可
        ### 注意！α和n_add是按sorted之后的顺序，而temp_two则是正常的顺序
        
        temp_two = cbind(temp_two,matrix(data = -1e12,ncol = max(n_add),nrow = k+1))
        for (i in 1:(k+1)){ #进行补充仿真,最后一行是上一代最优值
          if (n_add[i] != 0){
            if (i != k+1){
              for (j in (n_temp[i]+1):(n_temp[i]+n_add[i])){
                temp_two[pop_value_sorted[i,1],j] = object_func(pop[pop_value_sorted[i,1],1],pop[pop_value_sorted[i,1],2],pop[pop_value_sorted[i,1],3],pop[pop_value_sorted[i,1],4],pop[pop_value_sorted[i,1],5],pop[pop_value_sorted[i,1],6],pop[pop_value_sorted[i,1],7],pop[pop_value_sorted[i,1],8],pop[pop_value_sorted[i,1],9],pop[pop_value_sorted[i,1],10],sigma)
              }
            }else{
              for (j in (n_temp[i]+1):(n_temp[i]+n_add[i])){
                temp_two[i,j] = object_func(best_before[1],best_before[2],best_before[3],best_before[4],best_before[5],best_before[6],best_before[7],best_before[8],best_before[9],best_before[10],sigma)
              }
            }
          }else{
            next
          }
          if (i != k+1){
            pop_value[pop_value_sorted[i,1],2] = abs(mean(temp_two[pop_value_sorted[i,1],which(temp_two[pop_value_sorted[i,1],] != -1e12)]))
            pop_value[pop_value_sorted[i,1],3] = sd(temp_two[pop_value_sorted[i,1],which(temp_two[pop_value_sorted[i,1],] != -1e12)])
          }else{
            pop_value[i,2] = abs(mean(temp_two[i,which(temp_two[i,] != -1e12)]))
            pop_value[i,3] = sd(temp_two[i,which(temp_two[i,] != -1e12)])
          }
        }
        
        pop_value_sorted = rbind(pop_value[order(pop_value[-(k+1),2]),],pop_value[k+1,]) #新增计算量后按仿真值排序,最后一个仍然是上一代最优
        b = pop_value_sorted[1,1] #更新当代仿真值最优
        
        n_temp = n_temp+rep(max(n_add),k+1) #专用于储存仿真数据的矩阵列数
        
        ###注意，此时pop_value_sorted已经更新了排名，但是n和n_add的排列仍然是按旧的sorted顺序排列，因此需要更新n的顺序
        n[,2] = n[,2]+n_add
        temp_three = rep(NA,k+1)
        for (i in 1:k){
          temp_three[i] = n[which(n[,1] == pop_value_sorted[i,1]),2]
        }
        #View(temp_three)
        temp_three[k+1] = n[k+1,2]
        n = cbind(pop_value_sorted[,1],temp_three)
        
        #更新计算量
        sum_N = sum(n[,2])
        #View(alpha)
        
        #更新真实的α
        alpha = n[,2]/sum_N
        
        
        ############# t检验开始
        t_result = rep(NA,lim-1)
        df = rep(NA,lim-1)
        Q_set = matrix(data = NA, ncol=2,nrow = lim-1)
        for (i in 1:(lim-1)){
          t_result[i] = (pop_value_sorted[i+1,2]-pop_value_sorted[1,2])/sqrt(pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))
          df[i] = (pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/((pop_value_sorted[1,3]^2/(alpha[1]*sum_N))^2/(alpha[1]*sum_N-1)+(pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/(alpha[i+1]*sum_N-1))
          if (t_result[i] < qt(0.95,df = df[i])){
            Q_set[i,1] = pop_value_sorted[i+1,2]
            Q_set[i,2] = "Same"
          }else{
            Q_set[i,2] = "Not Same"
          }
        }
        
        Q_mean = sum(pop_value_sorted[c(1,which(Q_set[,2]=="Same")+1),2])/length(c(1,which(Q_set[,2]=="Same")+1)) ###更换均值
        ### 注意，这里的第一个是排名第二的个体与最优个体的对比
        
      }### 一代OCBA结束
      
      
      else{ #如果本代最优比上一代最优好 
        converge_rate = matrix(data = NA,nrow = k,ncol = 4) #用来保存三个值和类别
        for (i in 2:k){ #计算Ii+1,Ii-1,G
          if (i == 2){
            if (Q_set[i-1,2] == "Same"){
              converge_rate[i,1] = alpha[i]*(Q_mean-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(Q_mean+1e12)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (Q_mean-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            }else{
              converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]+1e12)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            }
          }else if (i == k){
            converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-1e12)^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
          }else if (i <= lim-1){
            if (Q_set[i-1,2] == "Same"){
              converge_rate[i,1] = alpha[i]*(Q_mean-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(Q_mean-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (Q_mean-pop_value_sorted[k+1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[k+1,3]^2/alpha[k+1]))
            }else{
              converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i+1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            }
          }else if (i == lim){
            if (Q_set[i-1,2] == "Same"){
              converge_rate[i,1] = alpha[i]*(Q_mean-M)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(Q_mean-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (Q_mean-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            }else{
              converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-M)^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,2] = alpha[i]*(pop_value_sorted[i,2]-pop_value_sorted[i-1,2])^2/(2*pop_value_sorted[i,3]^2)
              converge_rate[i,3] = (pop_value_sorted[i,2]-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            }
          }else{
            converge_rate[i,1] = alpha[i]*(pop_value_sorted[i,2]-M)^2/(2*pop_value_sorted[i,3]^2)
            converge_rate[i,2] = (pop_value_sorted[i,2]-pop_value_sorted[1,2])^2/(2*(pop_value_sorted[i,3]^2/alpha[i]+pop_value_sorted[1,3]^2/alpha[1]))
            converge_rate[i,3] = 1e12
          }
        }
        #View(alpha)
        
        sum_alpha = 0 #用于计算alpha[k+1]
        for (i in 2:k){ #先计算种群的alpha
          if (i <= lim-1){
            if (Q_set[i-1,2] == "Same"){
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i+1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }else{
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i+1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }
          }else if(i == lim){
            if (Q_set[i-1,2] == "Same"){
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-M))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(Q_mean-pop_value_sorted[1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }else{
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
                converge_rate[i,4] = 1
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-M))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
                converge_rate[i,4] = 2
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[i-1,2]))^2
              }
              if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 3){
                converge_rate[i,4] = 3
                alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[1,2]))^2
                sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
              }
            }
          }else{
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 1){
              converge_rate[i,4] = 1
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-M))^2
            }
            if (which(converge_rate[i,1:3] == min(converge_rate[i,1:3])) == 2){
              converge_rate[i,4] = 2
              alpha[i] = (pop_value_sorted[i,3]/(pop_value_sorted[i,2]-pop_value_sorted[1,2]))^2
              sum_alpha = sum_alpha+(alpha[i]/pop_value_sorted[i,3])^2
            }
          }
        }
        alpha[k+1] = (pop_value_sorted[k+1,3]/(pop_value_sorted[1,2]-pop_value_sorted[k+1,2]))^2
        alpha[1] = pop_value_sorted[1,3]*sqrt((alpha[k+1]/pop_value_sorted[k+1,3])^2+sum_alpha)
        
        #View(alpha)
        ### 开始增加计算量
        n_new = sum(c(n[,2],delta))/sum(alpha)*alpha
        n_add = rep(NA,k+1)
        for (i in 1:(k+1)){
          if (n_new[i] <= n[i,2]){
            n_add[i] = 0
          }
          else{
            n_add[i] = n_new[i]-n[i,2]
          }
        }
        n_add = round(n_add/sum(n_add)*delta)
        
        ### 新增计算量后的排名、μ、σ，得到当代最优b。此处pop_value_sorted的最后一个值一定是上一代最优
        ### 由于我们不需要更换上一代最优，因此直接增加仿真数量即可
        ### 注意！α和n_add是按sorted之后的顺序，而temp_two则是正常的顺序
        
        temp_two = cbind(temp_two,matrix(data = -1e12,ncol = max(n_add),nrow = k+1))
        for (i in 1:(k+1)){ #进行补充仿真,最后一行是上一代最优值
          if (n_add[i] != 0){
            if (i != k+1){
              for (j in (n_temp[i]+1):(n_temp[i]+n_add[i])){
                temp_two[pop_value_sorted[i,1],j] = object_func(pop[pop_value_sorted[i,1],1],pop[pop_value_sorted[i,1],2],pop[pop_value_sorted[i,1],3],pop[pop_value_sorted[i,1],4],pop[pop_value_sorted[i,1],5],pop[pop_value_sorted[i,1],6],pop[pop_value_sorted[i,1],7],pop[pop_value_sorted[i,1],8],pop[pop_value_sorted[i,1],9],pop[pop_value_sorted[i,1],10],sigma)
              }
            }else{
              for (j in (n_temp[i]+1):(n_temp[i]+n_add[i])){
                temp_two[i,j] = object_func(best_before[1],best_before[2],best_before[3],best_before[4],best_before[5],best_before[6],best_before[7],best_before[8],best_before[9],best_before[10],sigma)
              }
            }
          }else{
            next
          }
          if (i != k+1){
            pop_value[pop_value_sorted[i,1],2] = abs(mean(temp_two[pop_value_sorted[i,1],which(temp_two[pop_value_sorted[i,1],] != -1e12)]))
            pop_value[pop_value_sorted[i,1],3] = sd(temp_two[pop_value_sorted[i,1],which(temp_two[pop_value_sorted[i,1],] != -1e12)])
          }else{
            pop_value[i,2] = abs(mean(temp_two[i,which(temp_two[i,] != -1e12)]))
            pop_value[i,3] = sd(temp_two[i,which(temp_two[i,] != -1e12)])
          }
        }
        
        pop_value_sorted = rbind(pop_value[order(pop_value[-(k+1),2]),],pop_value[k+1,]) #新增计算量后按仿真值排序,最后一个仍然是上一代最优
        b = pop_value_sorted[1,1] #更新当代最优
        
        n_temp = n_temp+rep(max(n_add),k+1) #专用于储存仿真数据的矩阵列数
        
        ###注意，此时pop_value_sorted已经更新了排名，但是n的排列仍然是按旧的sorted顺序排列，因此需要更新n的顺序
        n[,2] = n[,2]+n_add
        temp_three = rep(NA,k+1)
        for (i in 1:k){
          temp_three[i] = n[which(n[,1] == pop_value_sorted[i,1]),2]
        }
        #View(temp_three)
        temp_three[k+1] = n[k+1,2]
        n = cbind(pop_value_sorted[,1],temp_three)
        
        #更新计算量
        sum_N = sum(n[,2])
        
        #更新真正的α
        alpha = n[,2]/sum_N
        
        
        ############# t检验开始
        t_result = rep(NA,lim-1)
        df = rep(NA,lim-1)
        Q_set = matrix(data = NA, ncol=2,nrow = lim-1)
        for (i in 1:(lim-1)){
          t_result[i] = (pop_value_sorted[i+1,2]-pop_value_sorted[1,2])/sqrt(pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))
          df[i] = (pop_value_sorted[1,3]^2/(alpha[1]*sum_N)+pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/((pop_value_sorted[1,3]^2/(alpha[1]*sum_N))^2/(alpha[1]*sum_N-1)+(pop_value_sorted[i+1,3]^2/(alpha[i+1]*sum_N))^2/(alpha[i+1]*sum_N-1))
          if (t_result[i] < qt(0.95,df = df[i])){
            Q_set[i,1] = pop_value_sorted[i+1,2]
            Q_set[i,2] = "Same"
          }else{
            Q_set[i,2] = "Not Same"
          }
        }
        
        Q_mean = sum(pop_value_sorted[c(1,which(Q_set[,2]=="Same")+1),2])/length(c(1,which(Q_set[,2]=="Same")+1)) ###更换均值
        ### 注意，这里的第一个是排名第二的个体与最优个体的对比
      }
    }
    ############# 至此OCBA结束，我们得到了混合种群的计算量分配n，当代最优b,排序后的种群信息pop_value_sorted(k+1行,3列。最后一行是上一代最优)
    
    
    ### 记录每代的最优个体,μ,真值
    history_best[z,] = c(pop[pop_value_sorted[1,1],],pop_value_sorted[1,2],pop_value_sorted[1,4])
    
    prob = rep(NA,k)
    for(i in 1:k){ #生成对应排名的权重
      prob[i] = (1/(q*k*sqrt(2*pi)))*exp(-(i-1)^2/(2*q^2*k^2))
    }
    cum_prob = cumsum(prob)/sum(prob)
    #View(cum_prob)
    
    ### Exploration阶段，随机选择m/2个体
    pop_chosen = matrix(data = NA,nrow = m,ncol = dimension+1)
    #View(cum_prob)
    for (i in 1:(m/2)){ #选择m个个体
      temp_four = runif(1)
      for (j in 1:k){
        if (j == 1 & temp_four <= cum_prob[j]){
          pop_chosen[i,] = c(pop_value_sorted[1,1],pop[pop_value_sorted[1,1],])
          break
        }else if(j == 1){
          next
        }else{
          if (temp_four <= cum_prob[j] & temp_four > cum_prob[j-1]){
            pop_chosen[i,] = c(pop_value_sorted[j,1],pop[pop_value_sorted[j,1],])
            break
          }
        }
      }
    }  
    
    ### 生成m/2个随机个体
    new_pop = matrix(data = NA,ncol = dimension+2,nrow = m) #前一半是exploration，后一半是exploitation;New pop还需要储存每个变量自身的方差
    for (i in 1:(m/2)){
      for (j in 1:(dimension)){
        stand_devia = zeta/(k-1)*sum(abs(pop[pop_chosen[i,1],j]-pop[,j]))
        new_pop[i,j] = rnorm(1,mean = pop_chosen[i,j+1],sd = stand_devia)
        new_pop[i,j+2] = stand_devia
        if (new_pop[i,j] > x_max){
          new_pop[i,j] = x_max
        }
        else if (new_pop[i,j] < x_min){
          new_pop[i,j] = x_min
        }
        else{
          next
        }
      }
    }
    
    
    ### Exploitation阶段，选择最优m/2个体并生成新个体
    for (i in 1:(m/2)){
      for (j in 1:dimension){
        stand_devia = zeta/(k-1)*sum(abs(pop[pop_value_sorted[i,1],j]-pop[,j]))
        new_pop[i+m/2,j] = rnorm(1,mean = pop[pop_value_sorted[i,1],j],sd = stand_devia)
        new_pop[i+m/2,j+2] = stand_devia
        if (new_pop[i+m/2,j] > x_max){
          new_pop[i+m/2,j] = x_max
        }
        else if (new_pop[i+m/2,j] < x_min){
          new_pop[i+m/2,j] = x_min
        }
        else{
          next
        }
      }
    }
    
    ### 精英引导阶段
    c = (z/max_iter)*(cf-ci)+ci
    for (i in 1:m){
      for (j in 1:dimension){
        if (pop_value_sorted[1,2] < pop_value_sorted[k+1,2]){ #如果当代最优更好，用当代最优引导
          if (abs(pop[b,j] - new_pop[i,j]) < 1e-06){
            new_pop[i,j] = new_pop[i,j] + c*runif(1,min = 0,max = 1)*new_pop[i,j+2]
            if (new_pop[i,j] > x_max){
              new_pop[i,j] = x_max
            }
            else if (new_pop[i,j] < x_min){
              new_pop[i,j] = x_min
            }else{
              next
            }
          }else{
            new_pop[i,j] = new_pop[i,j] + c*runif(1,min = 0,max = 1)*(pop[b,j]-new_pop[i,j])
            if (new_pop[i,j] > x_max){
              new_pop[i,j] = x_max
            }
            else if (new_pop[i,j] < x_min){
              new_pop[i,j] = x_min
            }else{
              next
            }
          }
        }else{ #如果当代最优不如上一代最优，用上一代最优引导
          if (abs(best_before[j] - new_pop[i,j]) < 1e-06){
            new_pop[i,j] = new_pop[i,j] + c*runif(1,min = 0,max = 1)*new_pop[i,j+2]
            if (new_pop[i,j] > x_max){
              new_pop[i,j] = x_max
            }
            else if (new_pop[i,j] < x_min){
              new_pop[i,j] = x_min
            }else{
              next
            }
          }else{
            new_pop[i,j] = new_pop[i,j] + c*runif(1,min = 0,max = 1)*(best_before[j]-new_pop[i,j])
            if (new_pop[i,j] > x_max){
              new_pop[i,j] = x_max
            }
            else if (new_pop[i,j] < x_min){
              new_pop[i,j] = x_min
            }else{
              next
            }
          }
        }
      }
    }
    
    ### 选择前k个种群，组成新种群
    hybird_pop = rbind(pop,new_pop[,1:dimension]) #将更新后的new_pop合并
    temp_five = matrix(data=NA,nrow=m,ncol=n1)
    new_pop_value = rep(NA,m)
    for (i in 1:m){ #计算new_pop未排序的均值
      for (j in 1:n1){
        temp_five[i,j] = object_func(new_pop[i,1],new_pop[i,2],new_pop[i,3],new_pop[i,4],new_pop[i,5],new_pop[i,6],new_pop[i,7],new_pop[i,8],new_pop[i,9],new_pop[i,10],sigma)
      }
      new_pop_value[i] = abs(mean(temp_five[i,]))
    }
    hybird_pop_value = cbind(seq(1,(m+k)),c(pop_value[-(k+1),2],new_pop_value)) #我们得到pop和new_pop的未排序的μ，其中new_pop的序号必定是k+1...k+m
    hybird_pop_value_sorted = hybird_pop_value[order(hybird_pop_value[,2]),] 
    for (i in 1:k){ #更新种群
      if (hybird_pop_value_sorted[i,1] > k){
        pop[i,] = new_pop[hybird_pop_value_sorted[i,1]-k,1:dimension]
      }else{
        pop[i,] = pop[hybird_pop_value_sorted[i,1],]
      }
    }
    
    
   # if (z >= 20){
  #   plot(seq(1,100),history_best[,4],type = "l")
   #  plot(seq(1,100),history_best[,3],type = "l")
   #  }
    
    
    
    #View(n)
    
    history_n[,z] = n[,2]
    z=z+1
    n = rep(n0,k+1)
  }
  #temp_six_avg[,a] = history_best[,3]
  #temp_six_true[,a] = history_best[,4]
  # }
  return(list(cbind(history_best[,(dimension+1):(dimension+2)]),history_n))
}

result = elite_new_partrank_ttest_OCBA(1)
View(result[[2]])


plot(seq(1,4000),result[[1]][,1],type = "l")





