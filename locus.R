
#################################
library(tidyverse)
library(glue)
require(data.table)
library(caret)
library(showtext)
library(ggokabeito)
#################################

font_add_google(name = "Caveat", family = "my_font")
showtext_auto()
theme_clean <- function() {
  theme_minimal(base_family = "my_font") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 36, family = "my_font"),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 40),
          axis.title = element_text(face = "bold", size =40),
          strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.text = element_text(size = 16))
}

resol <- 0.0003
dynamic_range <- 2 # 5
crowd_mag_limit <- 5
degree_decimals <- 4
M <- 0.1 # 0.1
dec.size = 10/60 # 0.167
dec.super = 10.5/60 # 0.172
'%+%' <- function(x,y) paste(x,y,sep="")

svr1 <- readRDS("models/svr-C20-sigma012-3_60mins")
preProcValues <- readRDS("models/preProcValues")

rev_logit <- function(x){
  x = x + preProcValues$mean["cor_logit"]
  x = x * preProcValues$std["cor_logit"]
  (exp(x)-1)/(exp(x)+1)
}

urlBase = "http://skyserver.sdss.org/dr15/SkyserverWS/SearchTools/SqlSearch?"



ObjID <- "1237680117417115655" # star
ObjID <- "1237674649391857845" #"1237648720693887182"
ObjID <- "1237651225171067249" # green publication QSO
# ObjID <- "1237651211213144306" # first QSO in publication 
ObjID <- "1237651224634851622" # second QSO in publication 
# SQL that downloads some info on the chosen target from SDSS.
# ObjID from SDSS specifies the target


locus_scoring <- function(index=1) {
  ObjID <- objid_list[index]
  
  targetSqlQuery <- paste("SELECT top 10 ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z FROM photoObj WHERE ObjID = ", ObjID) |> 
    str_squish() |> 
    str_replace_all(" ", "%20")
  
  
  # downloads target data
  # dataframe target has necessary info
  
  
  target <- read_csv(glue("{urlBase}cmd={targetSqlQuery}&format=csv"), skip = 1)
  
  # sets some variables for convenience. Last two are the field sizes
  # ra.size is automatically adjusted for each target depending on its dec
  # M is the maximum colour difference
  # resol is important for gauging crowded references
  # dynamic range is to prevent saturation of either target or reference
  u <- target$psfmag_u
  g <- target$psfmag_g
  r <- target$psfmag_r
  i <- target$psfmag_i
  z <- target$psfmag_z
  ra <- target$ra
  dec <- target$dec
  ra.size = dec.size / cos(dec*pi/180)
  ra.super = dec.super / cos(dec*pi/180)
  
  
  # SQL query that downloads data from SDSS for objects
  # potentially in the same field as the target
  mySqlQuery1 <-str_glue(
    "SELECT objID, ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z\n
  FROM photoObj\n
  WHERE (ra between ({ra - ra.size}) AND ({ra + ra.size})\n
  OR ra BETWEEN ({360 + ra - ra.size}) AND ({360 + ra + ra.size})\n
  OR ra BETWEEN ({-360 + ra - ra.size}) AND ({-360 + ra + ra.size})\n)
  AND dec BETWEEN ({dec - dec.size}) AND ({dec + dec.size})\n
  AND psfmag_r BETWEEN {r - dynamic_range} AND { r + dynamic_range} \n
  AND (psfmag_g - psfmag_r) BETWEEN ({g - r - M}) AND ({g - r + M})\n
  AND (psfmag_r - psfmag_i) BETWEEN ({r - i - M}) AND ({r - i + M})\n
  AND clean = 1") |> 
    str_squish() |> 
    str_replace_all(" ", "%20")
  
  # reads in data from SDSS.
  # dataframe called A has all the details
  
  A <- read_csv(glue("{urlBase}cmd={mySqlQuery1}&format=csv"), skip = 1) |> 
    mutate(objID = as.character(objID))
  
  # wrap-around for targets near 0 RA
  A$ra <- if_else(A$ra - target$ra > 180, A$ra - 360, A$ra) 
  # wrap-around for targets near 360 RA
  A$ra <- if_else(target$ra - A$ra > 180, A$ra + 360, A$ra)
  
  # function to calculate rating. Uses Oisin's routine
  # rating <- function(gr, rr, ir) {
  #   gt <- g
  #   rt <- r
  #   it <- i
  #   delta.CS <- (gt - rt) - (gr - rr)
  #   delta.CL <- (rt - it) - (rr - ir)
  #   RS <- 1 - abs(delta.CS / M)
  #   RL <- 1 - abs(delta.CL / M)
  #   RS * RL
  # }
  # 
  # 
  # # calculate ratings for each potential reference
  # ratings <- rating(A$psfmag_g, A$psfmag_r, A$psfmag_i)
  
  # ratings using Tom's SVR model
  tom_rating <- tibble(u = u, g = g, r = r, i = i, z = z,
                       u1 = A$psfmag_u,
                       g1 = A$psfmag_g,
                       r1 = A$psfmag_r,
                       i1 = A$psfmag_i,
                       z1 = A$psfmag_z,
                       cor_logit = 0.5)
  tom_transformed <- predict(preProcValues, tom_rating)
  tom_pred <- predict(svr1, tom_transformed) |> 
    rev_logit()
  ratings <- tom_pred^6
  
  # add ratings to the data frame
  A <- cbind(A, ratings)
  
  ##########################
  # finds all intersection points for each pair of potential references
  ##########################
  A_coords <- A %>% arrange(ra) %>% select(objID, ra, dec)
  
  swap_ij <- function(u, v){
    i <- v
    j <- u
  }
  
  int_pts_finder <- function(i, j) {
    if (i < j) {
      if (abs(A_coords$dec[i] - A_coords$dec[j]) < dec.size) {
        if (abs(A_coords$ra[i] - A_coords$ra[j]) < ra.size)
        {
          ifelse(
            A_coords$dec[i] > A_coords$dec[j],
            z <- data.frame(
              int_ra = c(A_coords$ra[i] + ra.size / 2,
                         A_coords$ra[j] - ra.size / 2),
              int_dec = c(A_coords$dec[j] + dec.size / 2,
                          A_coords$dec[i] - dec.size / 2),
              objID_i = c(A_coords$objID[i], A_coords$objID[i]),
              objID_j = c(A_coords$objID[j], A_coords$objID[j])
            ),
            z <- data.frame(
              int_ra = c(A_coords$ra[i] + ra.size / 2,
                         A_coords$ra[j] - ra.size / 2),
              int_dec = c(A_coords$dec[j] - dec.size / 2,
                          A_coords$dec[i] + dec.size / 2),
              objID_i = c(A_coords$objID[i], A_coords$objID[i]),
              objID_j = c(A_coords$objID[j], A_coords$objID[j])
            )
          )
          return(z)
        }
      }
    }
  }          
  index_matrix <- expand.grid(1:dim(A)[1], 1:dim(A)[1])
  names(index_matrix) <- c("i", "j")
  int.pts <- pmap(index_matrix, int_pts_finder) %>% 
    rbindlist() %>% 
    filter(between(int_ra, ra-ra.size/2, ra+ra.size/2), 
           between(int_dec, dec-dec.size/2, dec+dec.size/2)) %>%
    distinct()
  #####################################
  
  
  
  # function that returns the score for each intersection point
  score1 <- function(X, Y) {
    B <- A[abs(A$ra - X) <= ra.size / 2 + 0.001 & abs(A$dec - Y) <= dec.size / 2 + 0.001, ]
    sum(B$ratings)
  }
  
  # function that returns a dataframe with all references
  # in a FOV defined by an intersection point
  score2 <- function(X, Y) {
    A[abs(A$ra - as.numeric(X)) <= ra.size / 2 + 0.001 & abs(A$dec - as.numeric(Y)) <= dec.size / 2 + 0.001, ]
  }
  
  # calculates the score for each intersection point and orders them
  int.pts$score <- mapply(score1, int.pts$int_ra, int.pts$int_dec)
  int.pts <- int.pts[order(int.pts$score, decreasing = T), ]
  int.pts <- int.pts %>% filter(objID_i != ObjID, objID_j != ObjID)
  
  # prints out best pointing
  # usually a bunch of ties  but just picks the first one
  max.index <- which(int.pts$score == max(int.pts$score))
  
  # makes data frame of reference stars for best pointing
  B <- score2(int.pts[max.index[1], "int_ra"], int.pts[max.index[1], "int_dec"])
  
  C <- A |> 
    filter(between(ra, target$ra-ra.size/2, target$ra+ra.size/2)) |> 
    filter(between(dec, target$dec-dec.size/2, target$dec+dec.size/2))
  
  final_pointing <- data.frame(
    ObjID,
    ra = int.pts[max.index[1], ]$int_ra,
    dec = int.pts[max.index[1], ]$int_dec,
    score_locus = int.pts[max.index[1], ]$score,
    score_straight = sum(C$ratings),
    references_in_locus = nrow(B),
    references_in_straight = nrow(C),
    references_total = nrow(A)
  )
  
  write.csv(A, glue("objid-frames-all/objid{ObjID}.csv"))
  write.csv(B, glue("objid-frames-pointing/objid{ObjID}.csv"))
  
  final_pointing
}
locus_scoring_safely <- safely(locus_scoring)

my_stars <- read_csv("data/validate-525.csv") |>
  mutate(objid = as.character(objid),
         specobjid = as.character(specobjid),
         objid1 = as.character(objid1),
         specobjid1 = as.character(specobjid1))

objid_list <- my_stars$objid |> tail(325)

(t1 <- Sys.time())
z2 <- map(1:length(objid_list), ~locus_scoring_safely(.x)) %>%
  map_df("result") %>%
  compact()
(t2 <- Sys.time())
t2-t1
beepr::beep(5)
# 38s for first 20
# 13 minutes for first 200, 52 good

z2 |> 
  drop_na() |> 
  ggplot(aes(score_locus, score_straight)) + 
  geom_point(size = 5, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick4") +
  labs(x = "Score Using Locus Pointing",
       y = "Score Using Centered Target") +
  theme_clean()

write_csv(z2, "data/batch-of-tail-325.csv")
z3=read_csv("data/batch-of-head-200.csv") |> mutate(ObjID = as.character(ObjID))
z2=bind_rows(z2, z3)
