################# ԭʼOCBA������Ⱥ�㷨���� ##################
### �������
k=10 #kΪȺ��(design)����
dimension=2 #nΪ���ά��
zeta = 0.85 #��Ϊ������������
q = 0.2 #qΪȨ�ز���
m = 3 #mΪÿ�θ��µ�����
n0 = 10 #��ʼÿ��design�ļ�����
n = rep(n0,k) #�������design�ļ�����
N = 300 #������������ֵ
delta = 20 #ÿ�ε������Ӽ�������
sigma = 5 #��������С
object_func = function(x_one,x_two,sigma_test){ #Ŀ�꺯��
  #temp_one = (1+(x_one+x_two+1)^2*(19-14*x_one+3*x_one^2-14*x_two+6*x_one*x_two+3*x_two^2))*(30+(2*x_one-3*x_two)^2*(18-32*x_one+12*x_one^2+48*x_two-36*x_one*x_two+27*x_two^2))  Goldstein and Price �ɹ���
  #temp_one = x_one^2+2*x_two^2-(3/10)*cos(3*pi*x_one)-(2/5)*cos(4*pi*x_two)+7/10 B2 ������0��
  temp_one = (x_two-(5/(4*pi^2)*x_one^2)+5/pi*x_one-6)^2+10*(1-(1/(8*pi)))*cos(x_one)+10+rnorm(1,mean = 0,sd = sigma_test) #Branin
  return(temp_one)
}
object_func_true = function(x_one,x_two,sigma_test){ #Ŀ�꺯��
  temp_five = (x_two-(5/(4*pi^2)*x_one^2)+5/pi*x_one-6)^2+10*(1-(1/(8*pi)))*cos(x_one)+10 #Branin
  return(temp_five)
}
x_max = 15 #��������
x_min = -5 #��������

### ���ȷֲ����ɳ�ʼ����
pop = matrix(data = NA,nrow = k,ncol = dimension)
pop[,1] = runif(k,min = x_min,max = x_max)
pop[,2] = runif(k,min = x_min,max = x_max)

### ����ÿ��design�Ħ̺ͦҲ�����
pop_value = matrix(data = NA,nrow = k,ncol = 3) #3�зֱ𴢴���š��̡���
pop_value[,1] = seq(1,k)
temp_two = matrix(data = NA,ncol = n0,nrow = k)
for (i in 1:k){ #�����������
  for (j in 1:n0){
    temp_two[i,j] = object_func(pop[i,1],pop[i,2],sigma)
  }
}
for (i in 1:k){ #����̺ͦ�
  pop_value[i,2] = mean(temp_two[i,])
  pop_value[i,3] = sd(temp_two[i,])
}
pop_value_sorted = pop_value[order(pop_value[,2]),]
pop_value_sorted = cbind(pop_value_sorted,rep(NA,k))
for(i in 1:k){ #���ɶ�Ӧ������Ȩ��
  pop_value_sorted[i,4] = (1/(q*k*sqrt(2*pi)))*exp(-(i-1)^2/(2*q^2*k^2))
}

### ��¼��������壬����ֵ����ֵ
history_best = matrix(data = NA,nrow = 100,ncol = dimension+2)
history_best[1,] = c(pop[pop_value_sorted[1,1],],pop_value_sorted[1,2],object_func_true(pop[pop_value_sorted[1,1],1],pop[pop_value_sorted[1,1],2],sigma))


############# ACO������ʼ
ACO_OCBA = function(z_test){
  z = z_test
  while (z <= 100){
    ### ѡ�����
    pop_chosen = matrix(data = NA,nrow = m,ncol = dimension+1)
    cum_prob = cumsum(pop_value_sorted[,4])/sum(pop_value_sorted[,4])
    #View(cum_prob)
    for (i in 1:m){ #ѡ��m������
      temp_three = runif(1)
      for (j in 1:k){
        if (j == 1 & temp_three <= cum_prob[j]){
          pop_chosen[i,] = c(pop_value_sorted[1,1],pop[pop_value_sorted[1,1],])
          break
        }
        else if (j==1){
          next
        }
        else{
          if (temp_three <= cum_prob[j] & temp_three > cum_prob[j-1]){
            pop_chosen[i,] = c(pop_value_sorted[j,1],pop[pop_value_sorted[j,1],])
            break
          }
        }
      }
    }  
    
    ### ����m���������
    new_pop = matrix(data = NA,ncol = dimension,nrow = m)
    for (i in 1:m){
      for (j in 1:dimension){
        stand_devia = zeta/(k-1)*sum(abs(pop[pop_chosen[i,1],j]-pop[,j]))
        new_pop[i,j] = rnorm(1,mean = pop_chosen[i,j+1],sd = stand_devia)
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
    ### �������m��
    temp_four = pop_value_sorted[c(k,k-1,k-2),1]
    pop = rbind(pop[-temp_four,],new_pop) 
    
    ############# OCBA������ʼ
    ### ��������Ⱥ�ڼ�������Ϊn0�ĺ���ֵ���������õ���������b
    pop_value = matrix(data = NA,nrow = k,ncol = 3) #3�зֱ𴢴���š��̡���
    pop_value[,1] = seq(1,k)
    temp_two = matrix(data = NA,ncol = n0,nrow = k)
    for (i in 1:k){ #�����������
      for (j in 1:n0){
        temp_two[i,j] = object_func(pop[i,1],pop[i,2],sigma)
      }
      pop_value[i,2] = mean(temp_two[i,])
      pop_value[i,3] = sd(temp_two[i,])
    }
    pop_value_sorted = pop_value[order(pop_value[,2]),]
    b = pop_value_sorted[1,1] #��������
    sum_N = sum(n)
    
    ### ��ʼ�������Ӽ�����
    n_temp = n
    while (sum_N < N){
      ### ���������i
      ratio = rep(0,k)
      for (i in 1:k){ #������ķ������
        if (i!=b){
          ratio[i] = (pop_value[i,3]/(pop_value[b,2]-pop_value[i,2]))^2
        }
      }
      #View(ratio)
      ratio[b] = pop_value[b,3]*sqrt(sum(ratio/pop_value[,3])) #���Ž�ķ������
      
      ### ȷ������ļ�����
      n_new = rep(NA,k)
      for (i in 1:k){
        n_new[i] = sum(c(n,delta))*ratio[i]/sum(ratio)
      }
      n_add = rep(NA,k)
      for (i in 1:k){
        if (n_new[i] <= n[i]){
          n_add[i] = 0
        }
        else{
          n_add[i] = n_new[i]-n[i]
        }
      }
      n_add = round(n_add/sum(n_add)*delta)
      
      ### ������������ĺ���ֵ���������õ���������b
      temp_two = cbind(temp_two,matrix(data = -10000,ncol = max(n_add),nrow = k))
      for (i in 1:k){ #�����������
        if (n_add[i] != 0){
          for (j in (n_temp[i]+1):(n_temp[i]+n_add[i])){
            temp_two[i,j] = object_func(pop[i,1],pop[i,2],sigma)
          }
        }
        else{
          next
        }
        pop_value[i,2] = sum(temp_two[i,which(temp_two[i,] != -10000)])/length(temp_two[i,which(temp_two[i,] != -10000)])
        pop_value[i,3] = sd(temp_two[i,which(temp_two[i,] != -10000)])
      }
      pop_value_sorted = pop_value[order(pop_value[,2]),]
      b = pop_value_sorted[1,1] #������������ĵ�������
      
      n_temp = n_temp+rep(max(n_add),10) #ר���ڴ���������ݵľ�������
      n = n+n_add #���¼�����
      
      sum_N = sum(n)
    }
    ############# ���ˣ����ǵõ�������Ⱥ�ļ���������n����������b,��������Ⱥ��Ϣpop_value_sorted
    
    
    ### ���ɶ�ӦȨ��
    pop_value_sorted = cbind(pop_value_sorted,rep(NA,k))
    for(i in 1:k){
      pop_value_sorted[i,4] = (1/(q*k*sqrt(2*pi)))*exp(-(i-1)^2/(2*q^2*k^2))
    }
    
    ### ��¼ÿ�������Ÿ��������ֵ
    history_best[z,] = c(pop[b,],pop_value_sorted[1,2],object_func_true(pop[b,1],pop[b,2],sigma))
    
    z=z+1
  }
  plot(seq(1,100),history_best[,4],type = "l")
}

system.time(ACO_OCBA(2))









