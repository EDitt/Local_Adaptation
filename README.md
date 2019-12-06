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
- Number of Fruits (`Num_fruits`) - zero-truncated poisson distribution
  -The total number of fruits produced
  
Two other fitness components were analyzed in subsets of data:
- Number of seeds
  -
- Germination
  -Plants were grown from seed in 2014 and followed through reproduction
# Independent Variables
## Fixed Effects
- Soil type: Sandstone or Serpentine (2 levels) 
- Year: 2012-2015 (4 levels)
- Population: Sand or Serp (2 levels)
## Random Effects
- FlatEdge2: Whether or not a plant was transplanted into and edge or inner spot on a tray
- FieldEdge2: "" "" for field position (not necessarily the same as flat edge)
- Nested within year:
      NumDays_tray2field: The number of days between transplanting to tray vs field 
      Rep: a blocking factor. 2012 and 2013 only have 1 block per soil type. 2014 has       2 per soil type and 2015 has 4 per soil type

Will analyze 3 different fitness components-
1.) survival - modelled as bernoulli (0 or 1)
2.) fecundity (number of fruits) of individuals that survived to flower (Poisson)
3.) ASTER analyses that include survival, number of flowers, and number of fruits