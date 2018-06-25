################################################################
##############              MED 2018              ##############
##############           Adam Ma³kowski           ##############
##############   Implementacja algorytmu CLARA    ##############
################################################################

args = commandArgs(TRUE)

if(is.na(args[1]) | is.na(args[2])){
  cat('Brak koniecznych parametrów:\n
      [nazwa pliku]
      [liczba grup]
      [rozmiar próbki domyœlnie = 100]
      [liczba próbek domyœlnie = 10]
      [maks iteracji PAM domyœlnie = 5]
      [poziom metryki Minkowskiego domyœlnie = 2]')

} else{
  
  library(clusterCrit)
  library(dplyr)
  
  my_pam = function(data, n, max_iter = 5, minkowski_lvl = 2) {
    
    dist <- function(x1, x2) ((sum(abs(x1[ , !(names(x1) %in% c('cluster', 'u1', 'u2'))] - x2[ , !(names(x2) %in% c('cluster', 'u1', 'u2'))])) ^ minkowski_lvl))^(1/minkowski_lvl)
    
    calc_wjik = function(j, i, k) {
      ojk = o[as.numeric(row.names(data)[j]), 
              as.numeric(row.names(data)[k])]
      
      if(data$u1[j] < o[as.numeric(row.names(data)[j]), as.numeric(row.names(medoids)[i])]) {
        if(ojk >= data$u1[j]) {
          return(0)
        }
        else {
          return(ojk - data$u1[j])
        }
      }
      else {
        if(ojk < data$u2[j]) {
          return(ojk - data$u1[j])
        }
        else {
          return(data$u2[j]-data$u1[j])
        }
      }
    }
    
    update.data.frame = function() {
      for(i in 1:nrow(data)) {
        min_value = -1
        min_value_index = 0
        second_min_value = -1
        for(j in 1:nrow(medoids)) {
          dist = o[as.numeric(row.names(medoids)[j]), as.numeric(row.names(data)[i])]
          if(min_value < 0 | dist < min_value) {
            second_min_value = min_value
            min_value = dist
            min_value_index = j
          }
          else if(second_min_value < 0 | dist < second_min_value) {
            second_min_value = dist
          }
        }
        data$cluster[i] <<- min_value_index
        data$u1[i] <<- min_value
        data$u2[i] <<- second_min_value
      }
    }
    
    names = row.names(data)
    row.names(data) = seq(1,nrow(data),1)
    
    data$cluster = 0
    data$u1 = 0
    data$u2 = 0
    
    o = matrix(0, nrow = nrow(data), ncol = nrow(data))
    
    for(i in 1:nrow(data)) {
      for(j in 1:nrow(data)) {
        o[i, j] = dist(data[i,], data[j, ])
      }
    }
    medoids_indexes = sample(1:nrow(data), size=n, replace=FALSE)
    
    medoids = data[medoids_indexes, ]
    medoids$cluster = seq(1, n, 1)
    data = data[-medoids_indexes, ]
    
    update.data.frame()
    iter = 0
    while(iter < max_iter) {
      iter = iter + 1
      W = matrix(0, nrow = n, ncol = nrow(data), byrow = TRUE)
      
      for(k in 1:nrow(data)) {
        for(i in 1:n) {
          for(j in 1:nrow(data)) {
            W[i, k] = W[i, k] + calc_wjik(j, i, k)
          }
        }
      }
      minW = min(W)
      
      if(minW >= 0) {
        break
      }
      row = which(W==minW, arr.ind = TRUE)[1,][1]
      col = which(W==minW, arr.ind = TRUE)[1,][2]
      
      new_medoid = data[col, ] 
      old_medoid = medoids[row, ]
      
      medoids = rbind(medoids[-row, ], new_medoid)
      data = rbind(data[-col, ], old_medoid)
      
      medoids$cluster = seq(1, n, 1)
      
      update.data.frame()
    }
    medoids$clusters = NULL
    medoids$u1 = NULL
    medoids$u2 = NULL
    data$u1 = NULL
    data$u2 = NULL
    list(medoids, data$cluster)
  }
  
  my_clara = function(data, n, m = 100, k = 10, max_iter = 5, minkowski_lvl=2) {
    update.data.frame = function(data, medoids){
      sum = 0.0
      for(i in 1:nrow(data)) {
        min_value = -1
        min_value_index = 0
        for(j in 1:nrow(medoids)) {
          dist = dist(medoids[j, ], data[i, ])
          if(min_value < 0 | dist < min_value) {
            min_value = dist
            min_value_index = j
          }
        }
        sum = sum + min_value
        data$cluster[i] = min_value_index
      }
      quality <<- sum
      data
    }
    
    dist <- function(x1, x2) ((sum(abs(x1[ , !(names(x1) %in% c('cluster'))] - x2[ , !(names(x2) %in% c('cluster'))])) ^ minkowski_lvl))^(1/minkowski_lvl)
    
    data$cluster = 0
    min_quality = -1
    best_medoids = NULL
    for(i in 1:k) {
      cat('Clara - Iteracja:', i ,'\n')
      sample = sample_n(data[ ,!(names(data) %in% c('cluster'))], min(m, nrow(data)))
      medoids = my_pam(sample, n, max_iter, minkowski_lvl)[[1]]
      quality = 0
      result = update.data.frame(data, medoids)
      if(quality < min_quality | min_quality < 0) {
        cat('Clara - Znaleziono lepsze grupowanie!\n')
        min_quality = quality
        data = result
        best_medoids = medoids
      }
    }
    medoids$cluster = NULL
    list(medoids = medoids, clusters = data$cluster)
  }
  
  data = read.csv(args[1], header = FALSE)
  nums = sapply(data, is.numeric)

  n = as.integer(args[2])
  m = 100
  k = 10
  max_iter = 5
  minkowski_lvl = 2
  if(!is.na(args[3])){
    m = as.integer(args[3])
  }
  if(!is.na(args[4])){
    k = as.integer(args[4])
  }
  if(!is.na(args[5])){
    max_iter = as.integer(args[5])
  }
  if(!is.na(args[6])){
    minkowski_lvl = as.integer(args[6])
  }
  
  result = my_clara(data, n, m, k, max_iter, minkowski_lvl)
  indexes = intCriteria(as.matrix(mapply(data, FUN = as.numeric)), as.integer(result[[2]]), "all")
  
  capture.output(c(result, indexes), file = "clara_wyniki.txt")
  png('clara_wykres.png', width = 1000, height = 1000)
  plot(data, col = result[[2]])
  dev.off()
}
