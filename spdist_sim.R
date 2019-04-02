## A script to explore simulating species distributions using virtualspecies
#
# also using Tim's data to explore weighting pseudoabsences
#
# the goal is to test the weighting of pseudoabsences to a) improve land use
# models and b) improve edge-case uses of SDMS (e.g. fractured populations of
# a known-cosmopolitan species etc.)

# the first thing I want to try is using one of Tim's species that has a lot of
# data and cutting it up to reflect one of the ones that is very data poor.

# data here.

####### FUNCTIONS


make_bc_mask <- function(occurence, env_data) {
  # generate some background points - min 50.
  if (nrow(occurence) >= 50) {
    bg <- dismo::randomPoints(env, nrow(occurence))
  } else {
    bg <- dismo::randomPoints(env, 50)
  }

  occ_kf <- dismo::kfold(occurence, 5)
  bg_kf <- dismo::kfold(bg, 5)

  occ_train <- occurence[occ_kf != 1, ]
  occ_test <- occurence[occ_kf == 1, ]
  bg_train <- bg[bg_kf != 1, ]
  bg_test <- bg[bg_kf == 1, ]

  # fit model
  bc <- dismo::bioclim(env, occ_train)

  # evaluate
  bc_ev <- dismo::evaluate(occ_test, bg_test, bc, env)
  bc_t <- dismo::threshold(bc_ev, "spec_sens")

  # predict and mask
  x <- dismo::predict(env, bc)
  x <- x > bc_t
  ret <- list(model = bc, evaluation = bc_ev, mask = x)
  return(ret)
}

make_dat <- function(occurence, bg, env) {
  oc_kf <- kfold(occurence, 5)
  oc_train <- occurence[oc_kf != 1, ]
  oc_test <- occurence[oc_kf == 1, ]
  bg_kf <- dismo::kfold(bg, 5)
  bg_train <- bg[bg_kf != 1, ]
  bg_test <- bg[bg_kf == 1, ]
  train <- rbind(oc_train, bg_train)
  pa <- c(
    rep(1, nrow(oc_train)),
    rep(0, nrow(bg_train))
  )
  traindat <- raster::extract(env, train)
  traindat <- data.frame(cbind(pa = pa, traindat))
  test_oc <- data.frame(raster::extract(env, oc_test))
  test_bg <- data.frame(raster::extract(env, bg_test))

  ret <- list(
    traindat = traindat,
    test_pres = test_oc,
    test_bg = test_bg)
  return(ret)
}

sdm_both <- function(occurence, mask, env, model = "glm") {
  # Prepare MASKED absences
    mask[mask == 1] <- NA
    if (nrow(occurence) >= 50) {
      bg_m <- dismo::randomPoints(mask, nrow(occurence))
    } else {
      bg_m <- dismo::randomPoints(mask, 50)
    }

  # prepare UNMASKED absences
    if (nrow(occurence) >= 50) {
      bg_u <- dismo::randomPoints(env, nrow(occurence))
    } else {
      bg_u <- dismo::randomPoints(env, 50)
    }

  m <- stats::formula(paste("pa ~", paste(colnames(traindat_m)[-1], 
    collapse = "+")))
  dat_m <- make_dat(occurence, bg_m, env)
  dat_u <- make_dat(occurence, bg_u, env)

  if (model == "glm") {
    mod_m <- glm(formula(m), data = dat_m, 
      family = binomial(link = "logit"))
    mod_u <- glm(formula(m), data = dat_u, 
      family = binomial(link = "logit"))
  } else if (model == "rf") {
    mod_m <- suppressWarnings(randomForest::randomForest(m, data = dat_m,
      na.action = na.omit))
    mod_u <- suppressWarnings(randomForest::randomForest(m, data = dat_u,
      na.action = na.omit))
  }

  eval_m <- dismo::evaluate(test_oc, test_bg_m, mod_m)
  eval_u <- dismo::evaluate(test_oc, test_bg_u, mod_u)

  pred_m <- predict(env, mod_m)
  pred_u <- predict(env, mod_u)
  pred_m <- pred_m > threshold(eval_m, "spec_sens")
  pred_u <- pred_u > threshold(eval_u, "spec_sens")

  comparison <- c(masked_auc = eval_m@auc, unmasked_auc = eval_u@auc)

  ret <- list(
    comparison = comparison,
    masked = list(model = mod_m, evaluation = eval_m, prediction = pred_m),
    unmasked = list(model = mod_u, evaulation = eval_u, prediction = pred_u)
  )
  return(ret)
}


# setup the data to take a look...
library(raster)
library(sp)
library(dismo)

setwd("/home/hfg/Documents/hfg_soft/projects/wpa/species_and_layers")

spall <- read.csv("species_points.csv")

vnil <- vnil_sp <- spall[spall$spp == "Vachellia_nilotica", c("x", "y")]
sp::coordinates(vnil_sp) <- ~ x + y

# get the data etc.
env <- raster::stack(
  "./Env_layers/Vachellia_nilotica/bioclim_06_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_08_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_09_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_13_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_14_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_15_agg.asc",
  "./Env_layers/Vachellia_nilotica/bioclim_19_agg.asc",
  "./Env_layers/Vachellia_nilotica/lake_dist_resamp.asc",
  "./Env_layers/Vachellia_nilotica/slope_agg_10k.asc",
  "./Env_layers/Vachellia_nilotica/terrain_resamp.asc"
  )

names(env) <- c("bio06", "bio08", "bio09", "bio13", "bio14", "bio15", "bio19", 
  "lakes", "slope", "terrain")

# check the distribution of vnil
plot(env[[1]]);points(vnil)

# there are three obvious clumps - I think take two of those and cut out the
# rest - this will give the "widely distributed, very patchy" situation.

# take all the points south of -20, and then all the points west of 2? That
# should give the southerly points and the westerly points - those should be 
# two clumps...

vnil_s <- vnil[vnil$y <= -20, ]
vnil_w <- vnil[vnil$x <= 2, ]
vnil_train <- rbind(vnil_s, vnil_w)

vnil_e <- vnil[vnil$x > 2, ]
vnil_e <- vnil_e[vnil_e$y > -20, ]
vnil_test <- vnil_e

plot(env[[1]])
points(vnil_train, pch = 16, cex = 0.8, col = "black")
points(vnil_test, pch = 16, cex = 0.8, col = "red")

bc_mask <- make_bc_mask(vnil_train, env)

plot(bc_map_thresh)
points(vnil_train, pch = 16, cex = 0.8, col = "black")
points(vnil_test, pch = 16, cex = 0.8, col = "red")

# compare masked to unmasked performance
comp_glm <- sdm_both(vnil_test, bc_mask$mask, env, model = "glm")
comp_rf <- sdm_both(vnil_test, bc_mask$mask, env, model = "rf")

# plot predictions...
par(mfrow = c(1, 2))
plot(comp_glm$masked$prediction)
plot(comp_glm$unmasked$prediction)

par(mfrow = c(1, 2))
plot(comp_rf$masked$prediction)
plot(comp_rf$unmasked$prediction)

