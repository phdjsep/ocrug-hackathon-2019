library(rgdal)
library(PBSmapping)
library(MapGAM)

vio <- read.table("~/Desktop/Hackathon/CAnew/infractions_agg.csv", header = T, sep = ',')
vio[is.na(vio)] <- 0
View(vio)
library(dplyr)
vio12 <- filter(vio, issue_year == 2012)
vio13 <- filter(vio, issue_year == 2013)
vio14 <- filter(vio, issue_year == 2014)
vio15 <- filter(vio, issue_year == 2015)
vio16 <- filter(vio, issue_year == 2016)
vio17 <- filter(vio, issue_year == 2017)
vio18 <- filter(vio, issue_year == 2018)


CA <- read.table("~/Desktop/Hackathon/CAnew/CAnew.txt", header = T, sep = ',')
CA <- data.frame(CA)
CA <- CA[,c("NAME10", "X", "Y")]
CA$NAME10 <- toupper(CA$NAME10)
colnames(CA)[1] <- "county"

v12 <- merge(CA, vio12, by = "county")
v13 <- merge(CA, vio13, by = "county")
v14 <- merge(CA, vio14, by = "county")
v15 <- merge(CA, vio15, by = "county")
v16 <- merge(CA, vio16, by = "county")
v17 <- merge(CA, vio17, by = "county")
v18 <- merge(CA, vio18, by = "county")
View(v12)
View(v13)
View(v14)
View(v15)
View(v16)
View(v17)
View(v18)


CAmap <- readOGR(dsn = "~/Desktop/Hackathon/CAnew/CAnew.shp", layer="CAnew")
plot(CAmap)
points(CA$X, CA$Y, pch = 19)

CA.grid <- predgrid(X = CA$X, Y = CA$Y, map = CAmap)
CA.grid$pop <- 10000
dim(CA.grid)   # 5017   3
head(CA.grid)

plot(CAmap)
points(CA.grid, pch = 19)

gridlen <- nrow(CA.grid)
gridlen  # 5017

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v12)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v12)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v12[, 2:3]
dim(coords)
m.data <- v12
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.005

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v12.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v12.fit.c$fit)  # 3.616348e-05 1.308337e+02

colormap(v12.fit.c, CAmap, contours = "permrank", contours.lwd = 2, 
         legend.name = "Number of Violations per 10000 for 2012")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v12)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v12)
summary.glm(m.adj)
exp(-1.063e-04) # 0.9998937

anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v12[,2:3]
m.data <- v12
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.006

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v12.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v12.fit.a$fit)   # 1.610625e-01 9.550335e+04
range(v12.fit.c$fit)   # 8.030743e-06 2.195529e+02


colormap(v12.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 8.030743e-06, 
         mapmax = 2.195529e+02, legend.name = "Adjusted Number of Violations per 10000 in 2012" )

################################################### v13 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v13)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v13)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v13[, 2:3]
dim(coords)
m.data <- v13
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.005

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v13.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v13.fit.c$fit)  # 2.786487e-05 1.272587e+02

colormap(v13.fit.c, CAmap, contours = "permrank", contours.lwd = 2, 
         legend.name = "Number of Violations per 10000 for 2013")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v13)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v13)
summary.glm(m.adj)
exp(-1.063e-04) # 0.9998937

anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v13[,2:3]
m.data <- v13
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.002

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v13.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v13.fit.a$fit)   # 5.494455e-02 1.059628e+04
range(v13.fit.c$fit)   # 2.786487e-05 1.272587e+02


colormap(v13.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 5.494455e-02, 
         mapmax = 1.059628e+04, legend.name = "Adjusted Number of Violations per 10000 in 2013" )


################################################### v14 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v14)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v14)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v14[, 2:3]
dim(coords)
m.data <- v14
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.007

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v14.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v14.fit.c$fit) 

colormap(v14.fit.c, CAmap, contours = "permrank", contours.lwd = 2, 
         legend.name = "Number of Violations per 10000 for 2014")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v14)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v14)
summary.glm(m.adj)
exp(-1.063e-04) # 0.9998937

anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v14[,2:3]
m.data <- v14
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.002

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v14.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v14.fit.a$fit)   # 5.494455e-02 1.059628e+04
range(v14.fit.c$fit)   # 2.786487e-05 1.272587e+02


colormap(v14.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 5.494455e-02, 
         mapmax = 1.059628e+04, legend.name = "Adjusted Number of Violations per 10000 in 2014" )


################################################### v15 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v15)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v15)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v15[, 2:3]
dim(coords)
m.data <- v15
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.015

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v15.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v15.fit.c$fit)  # 2.477101e-04 8.327782e+01

colormap(v15.fit.c, CAmap, contours = "permrank", contours.lwd = 2, 
         legend.name = "Number of Violations per 10000 for 2015")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v15)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v15)
summary.glm(m.adj)


anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v15[,2:3]
m.data <- v15
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.001

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v15.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v15.fit.a$fit)   
range(v15.fit.c$fit)   # 2.477101e-04 8.327782e+01


colormap(v15.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 2.477101e-04, 
         mapmax = 8.327782e+01, legend.name = "Adjusted Number of Violations per 10000 in 2015" )



################################################### v16 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v16)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v16)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v16[, 2:3]
dim(coords)
m.data <- v16
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.009

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v16.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v16.fit.c$fit)  # 2.786487e-05 1.272587e+02

colormap(v16.fit.c, CAmap, contours = "permrank", contours.lwd = 2, legend.name = "Number of Violations per 10000 for 2016")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v16)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v16)
summary.glm(m.adj)


anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v15[,2:3]
m.data <- v15
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.001

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v16.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v16.fit.a$fit)   # 5.494455e-02 1.059628e+04
range(v16.fit.c$fit)   # 2.786487e-05 1.272587e+02


colormap(v16.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 5.494455e-02, 
         mapmax = 1.059628e+04, legend.name = "Adjusted Number of Violations per 10000 in 2016" )



################################################### v17 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v17)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v17)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v17[, 2:3]
dim(coords)
m.data <- v17
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.002

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v17.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v17.fit.c$fit)  # 7.094939e-06 2.813485e+01

colormap(v17.fit.c, CAmap, contours = "permrank", contours.lwd = 2, legend.name = "Number of Violations per 10000 for 2017")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v17)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v17)
summary.glm(m.adj)


anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v17[,2:3]
m.data <- v17
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.001

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v17.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v17.fit.a$fit)   # 2.121092e-03 4.264889e+02
range(v17.fit.c$fit)   # 7.094939e-06 2.813485e+01


colormap(v17.fit.a, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 7.094939e-06, 
         mapmax = 2.813485e+01, legend.name = "Adjusted Number of Violations per 10000 in 2017")



################################################### v18 ######################################################

m.crd <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = v18)
m.crd.0 <- gam(n_infractions_issued ~ offset(log(pop)), 
               family = poisson, data = v18)
anova(m.crd, m.crd.0)  # location matters

# Hypothesis Testing with Permutations
coords <- v18[, 2:3]
dim(coords)
m.data <- v18
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.crd.0, m.crd, test = "Chi")[2,4]
devrank[1]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
dim(permresults)   # 5017 1000

permresults[,1] <- predict(m.crd, CA.grid)

length(m.data$n_infractions_issued)
index <- sample(length(m.data$n_infractions_issued), replace = F)
m.data[,2:3] <- coords[index,]
m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
             family = poisson, data = m.data)
devrank[2] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]


i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[,2:3] <- coords[index, ]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)), 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.crd.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.crd <- (1000 - rank(devrank)[1])/1000
p.crd   # 0.016

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in 1:gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v18.fit.c <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.crd, CA.grid)), 
                  pointwise.permt = permrank[, 1]/1000)
range(v18.fit.c$fit)  # 3.677639e-03 3.887055e+06

colormap(v18.fit.c, CAmap, contours = "permrank", contours.lwd = 2, mapmin = 7.094939e-06, 
         mapmax = 2.813485e+01, legend.name = "Number of Violations per 10000 for 2018")



##################################### Adjusted Model ########################################

CA.grid$median_income <- 0

m.adj <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + median_income + offset(log(pop)), 
             family = poisson, data = v18)
m.adj.0 <- gam(n_infractions_issued ~ median_income + offset(log(pop)), 
               family = poisson, data = v18)
summary.glm(m.adj)


anova(m.adj.0, m.adj, test = "Chi")[2,5]  # 0

# Hypothesis Testing with Permutations
coords <- v18[,2:3]
m.data <- v18
devrank <- matrix(ncol = 1, nrow = 1000, 0)
devrank[1] <- anova(m.adj.0, m.adj, test = "Chi")[2, 4]
permresults <- matrix(ncol = 1000, nrow = gridlen, 0)
permresults[, 1] <- predict(m.adj, CA.grid)

i = j = 2
while (j < 1001){
  index <- sample(length(m.data$n_infractions_issued), replace = F)
  m.data[, 2:3] <- coords[index,]
  m.gam <- gam(n_infractions_issued ~ lo(X, Y, span = 0.15) + offset(log(pop)) + median_income, 
               family = poisson, data = m.data)
  devrank[i] <- anova(m.adj.0, m.gam, test = "Chi")[2,4]
  permresults[, i] <- predict(m.gam, CA.grid)
  if (j %% 100 == 0) print(i)
  i = j = j + 1
}

p.adj <- (1000 - rank(devrank)[1])/1000
p.adj   # 0.004

permrank <- matrix(ncol = 1000, nrow = gridlen, 0)

for (k in gridlen){
  permrank[k, ] <- rank(permresults[k, ])
}

v18.fit.a <- list(grid = CA.grid[, 1:2], fit = exp(predict(m.adj, CA.grid)), 
                 pointwise.permt = permrank[, 1]/1000)
range(v18.fit.a$fit)   # 3.677639e-03 3.887055e+06
range(v18.fit.c$fit)   # 3.677639e-03 3.887055e+06


colormap(v18.fit.a, CAmap, contours = "permrank", contours.lwd = 2,  legend.name = "Adjusted Number of Violations per 10000 in 2018")











