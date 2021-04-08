#### R script for Beck and Givnish, "Fine-scale environmental heterogeneity drives spatial niche partitioning among spring-flowering forest herbs"####

## R code by Jared Beck, Email: jared.j.beck@gmail.com, Website: jaredjbeck.github.io

## Load packages
library(ggplot2)
library(spatstat)
library(dplyr)
library(gridExtra)
library(gstat)
library(sp)
library(MuMIn)
library(emmeans)
library(nlme)
library(lme4)
library(car)

## Functions used in script

## Function to extract point-pattern simulation envelope for plotting
extract.envelope = function(env, from = NA, to = NA) {
  df = data.frame(distance = env$r,
                  observed = env$obs,
                  theoretical = env$theo,
                  lo = env$lo,
                  hi = env$hi,
                  from = from,
                  to = to)
  return(df)
}

## Function used to define breaks in spatial point pattern, used to visualize deviations from spatial randomness
define_breaks = function(xx, col = "col") {
  xx$is.break = NA
  xx$is.break[1] = 0
  for(i in 2:nrow(xx)){
    xx$is.break[i] = ifelse(xx[,col][i-1] == xx[,col][i], 0, 1)
  }
  
  breaks = which(xx$is.break == 1)
  
  xx$group2 = NA
  xx$group2[1:breaks[1]-1] = 1
  for(j in 1:(length(breaks)-1)){
    xx$group2[breaks[j]:breaks[j+1]] = j+1
  }
  xx$group2[breaks[length(breaks)]:nrow(xx)] = max(xx$group2, na.rm = T) +1
  return(xx)
}


#### Import data ####

## set working directory
setwd("C:/Users/jared/Documents/R_projects/natural_enemies_herbs/mcgilvra_mapping/for_submission/")

xx = read.csv("mcgilvra_herb_data.csv", stringsAsFactors = F)

soil = read.csv("mcgilvra_soil_data.csv", stringsAsFactors = F)

trees = read.csv("mcgilvra_tree_data.csv", stringsAsFactors = F)

## Tree community summary for focal plot
trees %>%
  group_by(species) %>%
  summarize(count = n(),
            basal.area.cm2 = sum(basal.area.cm2)) %>%
  group_by() %>%
  mutate(relative.basal.area = basal.area.cm2/sum(basal.area.cm2),
         basal.area.m2 = 0.0001* basal.area.cm2) %>%
  arrange(-relative.basal.area)

trees.small = trees %>%
  filter(basal.area.cm2 < pi*12.5^2) ## less than 25 cm DBH

trees.large = trees %>%
  filter(basal.area.cm2 >= pi*12.5^2) ## more than 25 cm DBH



## Calculate distances to nearest small/large trees
library(FNN)
xx$dist2smalltree = NA
xx$dist2largetree = NA


for(i in 1:dim(xx)[1]) {
  
  dist.df.small = rbind(xx[i,c("species","x", "y")],
                        trees.small[ ,c("species","x", "y")])
  
  xx$dist2smalltree[i] = knn.dist(data = dist.df.small[,c("x", "y")], k = 1, algorithm = "brute")[1]
  
  dist.df.large = rbind(xx[i,c("species","x", "y")],
                        trees.large[ ,c("species","x", "y")])
  
  xx$dist2largetree[i] = knn.dist(data = dist.df.large[,c("x", "y")], k = 1, algorithm = "brute")[1]
  
  
}## end loop

library(reshape2)

xx2 = xx %>%
  mutate(id = 1:length(site)) %>%
  melt(id.vars = c("id", "site", "species", "x", "y", "distance.to.tree.m", "soil.depth.cm")) %>%
  mutate(group = ifelse(variable %in% "dist2smalltree", "Small tree (<25 cm DBH)", "Large tree (>25 cm DBH"),
         distance = value)

ggplot(xx2, aes(x = species, y = distance, fill = group)) +
  geom_boxplot(position = position_dodge())



##Summary of sample sizes

xx %>% group_by(species) %>%
  summarize(N = n(),
            count_soil_50cm_plus = length(soil.depth.cm[soil.depth.cm > 49]),
            count_within_2m_of_tree = length(distance.to.tree.m[distance.to.tree.m < 2])) %>%
  mutate(percent.50plus = count_soil_50cm_plus/N,
         percent.2m = count_within_2m_of_tree/N)


ha = xx %>% filter(species %in% "Anemone acutiloba")
sc = xx %>% filter(species %in% "Sanguinaria canadensis")
tf = xx %>% filter(species %in% "Trillium flexipes")

#### << Figure S5. Relationship between soil depth and proximity to trees >> ####

with(xx, cor.test(soil.depth.cm, distance.to.tree.m))
with(xx, cor.test(soil.depth.cm, log(distance.to.tree.m)))

figs6 = ggplot(xx, aes(x = distance.to.tree.m, y = soil.depth.cm)) +
  geom_jitter(width = 0.1, height = 0, shape = 21, alpha = 0.5, fill = "grey50", color = "grey80") +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0,50,10), labels = c("0", "10", "20", "30", "40", "50+")) +
  stat_smooth(method = "lm", se = F, color = "black") +
  labs(x = "Distance to nearest tree (m)", y = "Soil depth (cm)") +
  annotate("text", x = 5, y = 3, hjust = 1, label = "r = 0.44, P < 0.001", fontface = "italic") +
  theme_classic()
figs6

# jpeg(filename = "mcgilvra_figS6.jpg", width = 4, height = 4, units = "in", res = 600, pointsize = 10)
# figs6
# dev.off()


#### << Figure 1. Map of herbs and trees >> ####
trees2 = trees
trees2$species = "trees"

trees2$species = factor(trees2$species, levels = c("Anemone acutiloba", "Sanguinaria canadensis", "Trillium flexipes", "trees"))
xx$species = factor(xx$species, levels = c("Anemone acutiloba", "Sanguinaria canadensis", "Trillium flexipes", "trees"))
trees2$basal.area.m2 = trees2$basal.area.cm2*0.0001

fig1 = ggplot(xx, aes(x= x, y = y, group = species, fill = soil.depth.cm, shape = species)) +
  coord_fixed() +
  facet_wrap(~species, ncol = 2) +
  geom_point(data = trees2, aes(x=x, y=y, size = basal.area.m2), shape = 1, color = "black", inherit.aes = F) +
  scale_size_continuous(name = expression(Basal~area~(m^2)), limits = c(0,1), breaks = c(0.1,0.25, 0.5, 0.75, 1)) +
  scale_shape_manual("Species", values = c(22, 23, 25)) +
  scale_x_continuous(limits = c(0,50))+ scale_y_continuous(limits = c(0,50)) +
  scale_fill_gradient("Soil depth (cm)", low = "white", high = "grey35", limits = c(-1,51), guide = guide_colorbar(reverse = T, frame.colour = "black", ticks.colour = "black")) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  geom_point(size = 1.5) +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 9, face = "italic"),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.justification = c(1,0),
        strip.text = element_text(face = "italic"))
fig1


pdf(file = "mcgilvra_fig1.pdf", width = 6, height = 5.25, pointsize = 10)
fig1
dev.off()



#### Comparing microsites occupied by each focal species ####

d = xx %>% select(-dist2smalltree, -dist2largetree )
d$soil.depth.cm = as.numeric(d$soil.depth.cm)

## Integrate background soil depth data with data for three focal species
dd = soil %>%
  mutate(site = "McGilvra Woods", 
         species = "background",
         soil.depth.cm = soil.depth.cm,
         distance.to.tree.m = NA) %>%
  select(site, species, x, y, distance.to.tree.m, soil.depth.cm)

dd = rbind(d, dd)

#### << Spatial ANOVAs for soil depth >> ####

## Semivariogram of soil depth
soil.vgm = variogram(soil.depth.cm~1, loc= ~x+y, data=dd, width = 0.25)
plot(soil.vgm)

## Exponential structure provides best fit
vgm.fit = fit.variogram(soil.vgm, vgm(c("Exp", "Sph", "Gau")))
vgm.fit
vgm.df = variogramLine(vgm.fit, maxdist = 25, n = 1000, min = 0.01)

fig3 = ggplot(vgm.df, aes(x = dist, y = gamma)) +
  geom_point(data = soil.vgm, aes(x = dist, y = gamma), size = 0.75, shape = 16, alpha = 0.5, inherit.aes = F) +
  geom_line() +
  scale_y_continuous("Semivariance of soil depth", limits = c(0,260), breaks = seq(0,250,50)) +
  scale_x_continuous("Distance between samples (m)") +
  annotate("text", x = 24, y = 1, hjust = 1, vjust = 0, label = "Nugget = 183.99\nPartial sill = 55.31\nRange = 1.72", size = 2.5) +
  theme_classic(base_size = 10)

pdf(file = "mcgilvra_fig3.pdf", width = 3.5, height = 3, pointsize = 10)
fig3
dev.off()


## Fitting ANOVAs and spatial ANOVAs

## Fit non-spatial Anova
sm = gls(soil.depth.cm ~ species, data = dd)

## Examine residual semivariance
plot(Variogram(sm, form = ~x + y, resType = "response", nint = 500, collapse = "fixed", maxDist = 25)) ## note spatial dependency within 2.5 m that is unaccounted for

## Fit models with different covariance structures
sm.exp <- gls( soil.depth.cm ~ species, correlation = corExp(form = ~x + y, nugget=T), data = dd )
sm.gauss <- gls( soil.depth.cm ~ species , correlation = corGaus(form = ~x + y, nugget=T), data = dd )
sm.spher <- gls( soil.depth.cm ~ species , correlation = corSpher(form = ~x + y, nugget=T), data = dd )

plot(Variogram(sm.spher, form = ~x + y, resType = "response", nint = 500, collapse = "fixed", maxDist = 25)) ## note the weak spatial signal in residuals
plot(Variogram(sm.gauss, form = ~x + y, resType = "response", nint = 500, collapse = "fixed", maxDist = 25)) ## note the weak spatial signal in residuals
plot(Variogram(sm.exp, form = ~x + y, resType = "response", nint = 500, collapse = "fixed", maxDist = 25)) ## note the weak spatial signal in residuals


AIC(sm, sm.exp, sm.gauss, sm.spher) ## the spherical correlation structure provides the best fit, though AIC for all 3 correlation structures within 2 AIC and all provide similar species-specific estimates

anova(sm.exp, sm) ## Including spatial covariance improves the model fit
anova(sm.exp)

## Create dataframe of estimated marginal means, bootstrapping takes a few minutes
semm.df = data.frame(emmeans(sm.exp, ~species,  mode = "df.error"))


#### << ANOVA for distance to nearest tree >> ####

tm = gls(distance.to.tree.m ~ species, data = d) ## non-spatial model
temm.df = data.frame(emmeans(tm, ~species))

anova(tm)


mod = lmer(distance ~ species*group + (1|id), data = xx2)
Anova(mod, test.statistic = "F")
emm.df = data.frame(emmeans(mod, ~species*group))


figS5 = ggplot(emm.df, aes(x = species, y  = emmean, shape = species, fill = group, group = group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.4), width = 0, size = 1) +
  geom_rect(aes(xmin = c(0.85,1.85,2.85, 1.05,2.05,3.05), xmax = c(0.95,1.95,2.95, 1.15,2.15,3.15), ymin = emmean - SE, ymax = emmean + SE), fill = "white", color = "black") +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  scale_shape_manual("Species", values = c(22, 23, 25, 22, 23, 25)) +
  scale_fill_grey() +
  annotate("text", x = Inf, y = -Inf, label = "Species: P < 0.001***\nTree size: P < 0.001***\nSpecies x Tree size: P = 0.002**",
           hjust = 1, vjust = -0.5, size = 2.5) +
  labs(x = "", y = "Estimated marginal mean\ndistance to nearest tree (m)") +
  scale_y_continuous(limits = c(0.5,4), breaks = seq(0,4,0.5)) +
  scale_x_discrete(labels = c("Anemone\nacutiloba", "Sanguinaria\ncanadensis", "Trillium\nflexipes")) +
  theme_classic(base_size = 9) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.position = "right",
        legend.text = element_text(face = "italic", size = 8),
        axis.text.x = element_text(face = "italic"),
        legend.title = element_blank())

jpeg(filename = "figS5.jpg", width = 5.5, height = 3.5, units = "in", res = 600, pointsize = 10)
figS5
dev.off()


#### << Figure 2. Comparison of microsites occupied by focal species >> ####

fig2a = ggplot(semm.df[!semm.df$species %in% "background",], aes(x = species, y  = emmean, shape = species)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = semm.df$lower.CL[semm.df$species %in% "background"],
  ymax = semm.df$upper.CL[semm.df$species %in% "background"]), alpha = 0.15, inherit.aes = F) +
  geom_hline(aes(yintercept = semm.df$emmean[semm.df$species %in% "background"]), linetype = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, size = 1) +
  geom_rect(aes(xmin = c(0.95,1.95,2.95), xmax = c(1.05,2.05,3.05), ymin = emmean - SE, ymax = emmean + SE), fill = "white", color = "black") +
  geom_point(size = 4, fill = "black") +
  annotate("text", x = 0.4, y = 42.7, label = "Soil depth of reference points", fontface = "italic", size = 3, hjust = -0.05, vjust = 0) +
  annotate("text", x = -Inf, y = Inf, label = "(a)", size = 4.5, hjust = -0.25, vjust = 1) +
  scale_shape_manual("Species", values = c(22, 23, 25)) +
  labs(x = "", y = "Estimated marginal mean\nsoil depth (cm)") +
  scale_x_discrete(labels = c("Anemone\nacutiloba", "Sanguinaria\ncanadensis", "Trillium\nflexipes")) +
  scale_y_continuous(limits = c(15,45), breaks = seq(15,45,10)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic"))

fig2b = ggplot(temm.df, aes(x = species, y  = emmean, shape = species)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, size = 1) +
  # geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), size = 2, color = "firebrick", width = 0) +
  geom_rect(aes(xmin = c(0.95,1.95,2.95), xmax = c(1.05,2.05,3.05), ymin = emmean - SE, ymax = emmean + SE), fill = "white", color = "black") +
  geom_point(size = 4, fill = "black") +
  annotate("text", x = -Inf, y = Inf, label = "(b)", size = 4.5, hjust = -0.25, vjust = 1) +
  scale_shape_manual("Species", values = c(22, 23, 25)) +
  labs(x = "", y = "Estimated marginal mean\ndistance to nearest tree (m)") +
  scale_y_continuous(limits = c(0.5,2.1), breaks = seq(0.5,2,0.5)) +
  scale_x_discrete(labels = c("Anemone\nacutiloba", "Sanguinaria\ncanadensis", "Trillium\nflexipes")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        # legend.position = c(0.95,0.05),
        legend.justification = c(1, 0),
        legend.background = element_rect(size = 0.5, color = "black"),
        legend.text = element_text(face = "italic", size = 8),
        axis.text.x = element_text(face = "italic"),
        legend.title = element_blank())

grid.arrange(fig2a, fig2b, ncol = 1)


pdf(file = "mcgilvra_fig2.pdf", width = 3.5, height = 6.5, pointsize = 10)
grid.arrange(fig2a, fig2b, ncol = 1)
dev.off()


#### Spatial point-pattern analysis using pair correlation functions ####


zz = xx
zz$group = zz$species
zz$distance.to.tree.m = NULL
zz$soil.depth.cm = NULL

trees$group = "tree"

## concatenate herb spatial data with tree coordinates
zz = rbind(zz[,c("site", "species", "x", "y", "group")], trees[,c("site", "species", "x", "y", "group")])

## create spatial objects for intraspecific point pattern analysis
ha.ppp = ppp(ha$x, ha$y, c(0,50), c(0,50))
tf.ppp = ppp(tf$x, tf$y, c(0,50), c(0,50))
sc.ppp = ppp(sc$x, sc$y, c(0,50), c(0,50))

## create spatial objects for interspecific/tree analyses
ha.tf.ppp = ppp(xx$x[xx$species %in% c("Anemone acutiloba", "Trillium flexipes")],
                xx$y[xx$species %in% c("Anemone acutiloba", "Trillium flexipes")],
                c(0,50), c(0,50))
marks(ha.tf.ppp) = xx[xx$species %in% c("Anemone acutiloba", "Trillium flexipes"), 'species']
ha.sc.ppp = ppp(xx$x[xx$species %in% c("Anemone acutiloba", "Sanguinaria canadensis")],
                xx$y[xx$species %in% c("Anemone acutiloba", "Sanguinaria canadensis")],
                c(0,50), c(0,50))
marks(ha.sc.ppp) = xx[xx$species %in% c("Anemone acutiloba", "Sanguinaria canadensis"), 'species']
sc.tf.ppp = ppp(xx$x[xx$species %in% c("Sanguinaria canadensis", "Trillium flexipes")],
                xx$y[xx$species %in% c("Sanguinaria canadensis", "Trillium flexipes")],
                c(0,50), c(0,50))
marks(sc.tf.ppp) = xx[xx$species %in% c("Sanguinaria canadensis", "Trillium flexipes"), 'species']

ha.tree.ppp = ppp(zz$x[zz$group %in% c("Anemone acutiloba", "tree")],
                  zz$y[zz$group %in% c("Anemone acutiloba", "tree")],
                  c(0,50), c(0,50))
marks(ha.tree.ppp) = zz[zz$group %in% c("Anemone acutiloba", "tree"), 'group']

sc.tree.ppp = ppp(zz$x[zz$group %in% c("Sanguinaria canadensis", "tree")],
                  zz$y[zz$group %in% c("Sanguinaria canadensis", "tree")],
                  c(0,50), c(0,50))
marks(sc.tree.ppp) = zz[zz$group %in% c("Sanguinaria canadensis", "tree"), 'group']

tf.tree.ppp = ppp(zz$x[zz$group %in% c("Trillium flexipes", "tree")],
                  zz$y[zz$group %in% c("Trillium flexipes", "tree")],
                  c(0,50), c(0,50))
marks(tf.tree.ppp) = zz[zz$group %in% c("Trillium flexipes", "tree"), 'group']


## Estimate inhomogeneous pair correlation functions and simulation envelopes
n.sim = 499
ha2ha = envelope.ppp(ha.ppp, pcfinhom, nsim = n.sim, correction = "Ripley")
sc2sc = envelope.ppp(sc.ppp, pcfinhom, nsim = n.sim, correction = "Ripley")
tf2tf = envelope.ppp(tf.ppp, pcfinhom, nsim = n.sim, correction = "Ripley")

ha2tf = envelope.ppp(ha.tf.ppp, pcfcross.inhom, nsim = n.sim, i = "Anemone acutiloba", j = "Trillium flexipes", correction = "Ripley")
tf2ha = envelope.ppp(ha.tf.ppp, pcfcross.inhom, nsim = n.sim, i = "Trillium flexipes", j = "Anemone acutiloba", correction = "Ripley")
ha2sc = envelope.ppp(ha.sc.ppp, pcfcross.inhom, nsim = n.sim, i = "Anemone acutiloba", j = "Sanguinaria canadensis", correction = "Ripley")
sc2ha = envelope.ppp(ha.sc.ppp, pcfcross.inhom, nsim = n.sim, i = "Sanguinaria canadensis", j = "Anemone acutiloba", correction = "Ripley")
sc2tf = envelope.ppp(sc.tf.ppp, pcfcross.inhom, nsim = n.sim, i = "Sanguinaria canadensis", j = "Trillium flexipes", correction = "Ripley")
tf2sc = envelope.ppp(sc.tf.ppp, pcfcross.inhom, nsim = n.sim, i = "Trillium flexipes", j = "Sanguinaria canadensis", correction = "Ripley")

ha2tree = envelope.ppp(ha.tree.ppp, pcfcross.inhom, nsim = n.sim, i = "Anemone acutiloba", j = "tree", correction = "Ripley")
sc2tree = envelope.ppp(sc.tree.ppp, pcfcross.inhom, nsim = n.sim, i = "Sanguinaria canadensis", j = "tree", correction = "Ripley")
tf2tree = envelope.ppp(tf.tree.ppp, pcfcross.inhom, nsim = n.sim, i = "Trillium flexipes", j = "tree", correction = "Ripley")


#### << Figure 3. Visualizing spatial patterns from inhomogeneous poisson process >> ####

## extract simulation envelopes and turn into data frame for plotting
pcf.i.df = rbind(extract.envelope(ha2ha, from = "Anemone", to = "Anemone"),
                 extract.envelope(sc2sc, from = "Sanguinaria", to = "Sanguinaria"),
                 extract.envelope(tf2tf, from = "Trillium", to = "Trillium"),
                 extract.envelope(ha2tf, from = "Anemone", to = "Trillium"),
                 extract.envelope(ha2sc, from = "Anemone", to = "Sanguinaria"),
                 extract.envelope(sc2tf, from = "Sanguinaria", to = "Trillium"),
                 extract.envelope(ha2tree, from = "Anemone", to = "trees"),
                 extract.envelope(sc2tree, from = "Sanguinaria", to = "trees"),
                 extract.envelope(tf2tree, from = "Trillium", to = "trees")) %>%
  filter(distance >= 0.1) %>% ## Exclude zero distance
  mutate(pattern = ifelse(observed > hi, "Aggregated",
                          ifelse(observed < lo, "Segregated", "Random"))) %>%
  mutate(group = paste(from, "\U2192", to, sep = " ")) %>%
  group_by(group) %>%
  mutate(max = max(c(observed, hi)))

pcf.i.df = data.frame(pcf.i.df)
pcf.i.df$pattern = factor(pcf.i.df$pattern, levels = c("Aggregated", "Random", "Segregated"))

## create bar to more easily display deviations from spatial randomness
pcf.i.bar = define_breaks(pcf.i.df, col = "pattern") %>%
  group_by(group) %>%
  mutate(ymax = max(c(observed, hi))) %>%
  group_by(group2, pattern, group, ymax) %>%
  summarize(xmin = min(distance), xmax = max(distance))


pcf.i.df$group = factor(pcf.i.df$group)
pcf.i.df$group = factor(pcf.i.df$group, 
       levels = c("Anemone → Anemone", "Sanguinaria → Sanguinaria", "Trillium → Trillium",
                  "Anemone → Sanguinaria", "Anemone → Trillium", "Sanguinaria → Trillium",
                  "Sanguinaria → Anemone", "Trillium → Anemone", "Trillium → Sanguinaria",
                  "Anemone → trees", "Sanguinaria → trees", "Trillium → trees"))

pcf.label.df = data.frame(group = c("Anemone → Anemone", "Sanguinaria → Sanguinaria", "Trillium → Trillium",
                                    "Anemone → Sanguinaria", "Anemone → Trillium", "Sanguinaria → Trillium",
                                    "Anemone → trees", "Sanguinaria → trees", "Trillium → trees"),
                          lab = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"))
pcf.label.df$group = factor(pcf.label.df$group, 
                        levels = c("Anemone → Anemone", "Sanguinaria → Sanguinaria", "Trillium → Trillium",
                                   "Anemone → Sanguinaria", "Anemone → Trillium", "Sanguinaria → Trillium",
                                   "Sanguinaria → Anemone", "Trillium → Anemone", "Trillium → Sanguinaria",
                                   "Anemone → trees", "Sanguinaria → trees", "Trillium → trees"))


pcf.i.bar$group = factor(pcf.i.bar$group)
pcf.i.bar$group = factor(pcf.i.bar$group, 
                        levels = c("Anemone → Anemone", "Sanguinaria → Sanguinaria", "Trillium → Trillium",
                                   "Anemone → Sanguinaria", "Anemone → Trillium", "Sanguinaria → Trillium",
                                   "Sanguinaria → Anemone", "Trillium → Anemone", "Trillium → Sanguinaria",
                                   "Anemone → trees", "Sanguinaria → trees", "Trillium → trees"))

#### << Figure 4. Pair correlation functions >> ####
fig4 = ggplot(pcf.i.df, aes(x = distance, y = observed)) +
  facet_wrap(~group, scales = "free", ncol = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(color = "black", size = 0.5) +
  labs(x = "Lag distance (m)", y = expression(italic(g[inhom]~(r)))) +
  geom_rect(data = pcf.i.bar, aes(xmin = xmin, xmax = xmax, ymin = 0-ymax/15, ymax = 0-ymax/60, color = pattern, fill = pattern), inherit.aes = F) +
  geom_rect(aes(xmin = 0, xmax = 12.5, ymin = 0-max/15, ymax = 0-max/60), color = "grey10", fill = NA, size = 0.5) +
  scale_fill_brewer(type = "div", palette = 7, name = "Spatial pattern") +
  scale_color_brewer(type = "div", palette = 7, name = "Spatial pattern") +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(dat = pcf.label.df, aes(x = 0.5, y = Inf, label = lab), hjust = 0, vjust = 1.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(limits = c(-0.1, 12.6), breaks = seq(0,12.5, 2.5)) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99,0.99),
        legend.justification = c(1,1),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.75,"line"),
        legend.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(face = "italic"))
print(fig4)

pdf(file = "mcgilvra_fig4.pdf", width = 7.25, height = 6, pointsize = 10)
fig4
dev.off()
