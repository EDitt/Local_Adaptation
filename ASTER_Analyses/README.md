# Local_Adaptation
Analyses for Local Adaptation manuscript

---

Reciprocal transplant experiments were performed in both parental habitats (Serpentine and Sandstone) from 2012-2015.
# Fitness Components
There were three main fitness components included in ASTER analyses:
- Survival (`Surv_flr`) - modelled as bernoulli (0 or 1)
  -whether or not the plant survived to flower
- Number of Flowers (`Num_flrs`) - zero-truncated poisson distribution
  -The total number of flowers produced
- Number of Fruits (`Num_fruits`) - poisson distribution
  -The total number of fruits produced

# Independent Variables
## Fixed Effects
- Soil type: Sandstone or Serpentine (2 levels) 
- Year: 2012-2015 (4 levels)
- Source Population: Sand or Serp (2 levels)  
- Edge: Whether or not a plant was transplanted to the outer edge of a field plot (2 levels: edge vs. non-edge)

## Random Effects
- In 2014-2015 there were multiple plots per soil type ("Rep"). These were included as random effects only in analyses of these years individual (including them in the full analyses obscured variation due to the year effect because there was only 1 rep per soil type in 2012-2013)
 
# Other fitness components
Two other fitness components were analyzed in subsets of data:
- Number of seeds
  -Mature fruits were collected from plants over the four years of data collected and seeds per fruit counted
- Germination
  -Plants were grown from seed in 2014 and followed through reproduction
