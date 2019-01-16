dgp= function(i,n, intercept.var, slope.var,cor, eb_0, eb_1, eb_2, sill, range){
  
  ###Set up an exponential covariance model:######
  sim <-  grf(n, cov.pars=c(sill, range), cov.model="exponential")
  
  #extract coordinates and convert into a spatialpointsdataframe object:
  coords<- data.frame(sim$coords[,1], sim$coords[,2])
  sim.spdf= SpatialPointsDataFrame(coords, data=data.frame(sim$data))
  
  #create vectors of the x and y points with the following spacing:
  x <- seq(bbox(sim.spdf)[1,1], bbox(sim.spdf)[1,2], length.out=sqrt(i))
  y <- seq(bbox(sim.spdf)[2,1], bbox(sim.spdf)[2,2], length.out=sqrt(i))
  
  #Creating a grid of all the pairs of coordinates and converting to a gridded object:
  xy <- expand.grid(x = x, y = y)
  grid.pts<-SpatialPointsDataFrame(coords= xy, data=xy)
  gridded(grid.pts) <- TRUE
  gridded(grid.pts)
  
  #Convert the grid of points into a SpatialPolygons:
  grid <- as(grid.pts, "SpatialPolygons")
  
  #Converting into a SpatialPolygonsDataFrame:
  gridspdf <- SpatialPolygonsDataFrame(grid, data=data.frame(block=row.names(grid), row.names=row.names(grid)))
  names.grd<-sapply(gridspdf@polygons, function(x) slot(x,"ID"))
  
  #Using over to identify which point belongs to which block:
  overlaid.pts = over(sim.spdf,gridspdf)
  
  block.pts = cbind(sim.spdf@data, overlaid.pts)
  block.pts$block= gsub("g", "", as.character(block.pts$block))
  block.pts$block= as.factor(block.pts$block)
  
  
  ####Simulating multilevel data####
  
  #i blocks/groups
  
  df <- data.frame(block = c(1:i))
  
  #true random intercept variance= 1, true random slope variance= 4
  cov.matrix <-  matrix(c(intercept.var^2, intercept.var*slope.var*cor, intercept.var*slope.var*cor, slope.var^2), nrow = 2,   byrow = TRUE)    
  
  #block-level random effects for intercept and slope:
  random.effects <-  rmvnorm(i, mean = c(0, 0), sigma = cov.matrix)
  
  #Adding noise and constructing the  reg coefficient for intercept (block-level):
  df$b0 = eb_0  +  random.effects[, 1]                #random intercept
  df$b1 = eb_1 + random.effects[, 2]                #individual-level covariate with a random slope
  df$b2=  eb_2                                           #group-level covariate
  
  #Generating group-level covariate "x2" such that each point belonging to the same block has the same group-level covariate value:  
  within_block= block.pts %>% group_by(block) %>%  mutate(individual = sequence(n()))
  within_block=within_block %>% group_by(block) %>% mutate(x2=rnorm(1))
  
  #Merging the df dataset (group-level info) and the within_block dataset:
  df = merge(df, within_block)
  
  #Generating an individual-level covariate "x1":
  df=df %>% mutate(x1= rnorm(n), e=rnorm(n))
  
  #Creating Y: 
  df$Y1 <- df$b0 + df$x1*df$b1 + df$x2*df$b2 + df$e
  
  #Running k-nearest neighbours for k=4, and testing for spatial autocorrelation with Moran's I (we expect no   autocorrelation since no spatial dependence has been incorporated yet)
  neighbours=knearneigh(as.matrix(sim$coords), 4)
  moran1= moran.mc(df$Y1, listw=nb2listw(knn2nb(neighbours)), 1000)
  
  #Adding in spatial dependence and testing for spatial autocorrelation:
  df$Y2 <- df$b0 + df$x1*df$b1 + df$x2*df$b2 + df$e + df$sim.data
  moran2= moran.mc(df$Y2, listw=nb2listw(knn2nb(neighbours)), 1000)
  
  
  #Running a multilevel model with a random intercept and a random slope on Y1:
  m0= (lmer(df$Y1~ x1 + x2 + (1|block) + (0+ x1|block), data=df))
  
  
  #Running a multilevel model with a random intercept and a random slope:
  m1= (lmer(df$Y2~ x1 + x2 + (1|block) + (0+ x1|block), data=df))
  
  #Running a CAR model:
  s1= spautolm(df$Y2~ x1 + x2, data = df, nb2listw(knn2nb(neighbours)),  family = "CAR")
  
  
  #Estimated coefficients for m0:
  estimated_coeff_m0= c(coef(summary(m0))["(Intercept)" , "Estimate"],
                        coef(summary(m0))["x1" , "Estimate"],
                        coef(summary(m0))["x2" , "Estimate"])
  
  
  #Estimated coefficients for m1:
  estimated_coeff_m1= c(coef(summary(m1))["(Intercept)" , "Estimate"],
                        coef(summary(m1))["x1" , "Estimate"],
                        coef(summary(m1))["x2" , "Estimate"])
  
  #Estimated coefficients for s1:
  estimated_coeff_s1= c(coef(summary(s1))["(Intercept)"],
                        coef(summary(s1))["x1"],
                        coef(summary(s1))["x2"])
  
  #True coefficients: 
  true_coeff=c(eb_0,  eb_1,  eb_2)
  
  #table to compare coefficients
  coeff_df= data.frame(estimated_coeff_m0, estimated_coeff_m1, estimated_coeff_s1, true_coeff, row.names = c("Intercept", "x1", "x2"))  
  
  #Estimated variance components for m0:
  estimated_var_m0= c(data.frame(VarCorr(m0),comp="Variance")[1,4],   #intercept
                      data.frame(VarCorr(m0),comp="Variance")[2,4],  #x2  
                      data.frame(VarCorr(m0),comp="Variance")[3,4])   #residual
  
  #Estimated variance components for m1:
  estimated_var_m1= c(data.frame(VarCorr(m1),comp="Variance")[1,4],   #intercept
                      data.frame(VarCorr(m1),comp="Variance")[2,4],  #x2  
                      data.frame(VarCorr(m1),comp="Variance")[3,4])   #residual
  
  #Estimated residual variance for s1:
  estimated_var_s1= c(NA, NA, (s1[["fit"]][["s2"]]))   #residual
  
  #true variance:
  true_var= c(intercept.var^2, slope.var^2,1)
  
  #table to compare variance
  var_df= data.frame(estimated_var_m0, estimated_var_m1, estimated_var_s1, true_var, row.names = c("Intercept", "x2", "residual"))
  
  #Moran df:
  moran_df = data.frame(moran1[["p.value"]], moran2[["p.value"]], row.names= c("p.value"))
  
  return(list(points(sim) ,plot(gridspdf, add=T), coeff_df, var_df, moran_df))}




  
  