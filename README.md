# Demographic models predicting the effect of human hunting on now-extinct dwarf hippopotamus 🦛 and dwarf elephant 🐘 in Cyprus 🇨🇾
<a href="https://www.ucy.ac.cy/migrate/"><img align="right" src="www/MIGRATElogo.jpg" width="200" style="margin-top: 20px"></a>

Part of the <a href="https://www.ucy.ac.cy/migrate/">MIGRATE</a> (<strong>M</strong>odell<strong>i</strong>ng Demo<strong>gr</strong>aphy and <strong>A</strong>daptation in the Initial Peopling of the Eastern M<strong>e</strong>diterranean Islandscape) project, under the auspices of the European Union Research and Innovation Foundation for Research, Technological Development and Innovation "Restart 2016-2020".
<br>
<br>
<strong>lead investigator</strong>: Dr <a href="https://ucy.academia.edu/TheodoraMoutsiou">Theodora Moutsiou</a><br>
<strong>key personnel</strong>: Dr <a href="https://scholar.google.com.au/citations?user=BU25ogMAAAAJ&hl=en">Christian Reepmeyer</a>, Associate Professor <a href="https://www.ucy.ac.cy/directory/en/profile/demest">Stella Demesticha</a>, Dr <a href="https://www.ucy.ac.cy/directory/en/profile/arkasian">Vasiliki Kassianidou</a>, Dr <a href="https://www.cut.ac.cy/faculties/fet/ceg/staff/athos.agapiou/?languageId=1">Athos Agapiou</a>, Dr <a href="https://www.researchgate.net/profile/Zomenia-Zomeni">Zomenia Zomeni</a>, Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
<strong>collaborators</strong>: Dr <a href="https://globalecologyflinders.com/people/#COORDINATOR">Frédérik Saltré</a>, Dr <a href="https://scholar.google.com.au/citations?user=-BSGg1MAAAAJ&hl=en">Salvador Herrando-Pérez</a>
<br>
## Project overview
Project <a href="https://www.ucy.ac.cy/migrate/">MIGRATE</a> seeks to offer novel insights into population dynamics and range shifts that resulted in dispersals from the Eastern Mediterranean mainland to the island of <a href="https://www.google.com/maps/place/Cyprus/@35.1670135,32.765821,9z/">Cyprus</a> at a critical period (Late Pleistocene, 45-12 ka) through stochastic spatial modelling. This advanced modelling will  enhance our understanding of timing, and climatic and social factors important in the initial colonisation of Cyprus. The proposed project aims to establish new research domains in the field of Cypriot archaeology extending traditional chronological frontiers beyond the Holocene (current Warm Period), encompassing innovative and interdisciplinary methodologies at the forefront of archaeological research.

## Papers arising
- Bradshaw, CJA, C Reepmeyer, F Saltré, A Agapiou, V Kassinadiou, S Demesticha, Z Zomeni, M Polidorou, T Moutsiou. 2023. <a href="http://doi.org/10.21203/rs.3.rs-3468157/v1">Demographic models predict end-Pleistocene arrival and rapid expansion of pre-agropastoralist humans in Cyprus</a>. <em><strong>Research Square</strong></em> (pre-print) doi:10.21203/rs.3.rs-3468157/v1 (see also <a href="https://github.com/cjabradshaw/CyprusHumanPleistocene">related Github repository</a>)
- Moutsiou T. 2021. <a href="http://doi.org/10.1016/j.quaint.2020.09.012">Climate, environment and cognition in the colonisation of the Eastern Mediterranean islands during the Pleistocene</a>.  <em><strong>Quaternary International</strong></em> 577:1-14
- Moutsiou T, C Reepmeyer, V Kassianidou, Z Zomeni, A Agapiou A. 2021. <a href="http://doi.org/10.1371/journal.pone.0258370">Modelling the Pleistocene colonisation of Eastern Mediterranean islandscapes</a>. <em><strong>PLoS One</strong></em> 16:e0258370

## Focal manuscript
Bradshaw, CJA, F Saltré, S Herrando-Pérez, C Reepmeyer, T Moutsiou. Palaeolithic human populations on Cyprus and the mechanisms of hunting native megafauna to extinction (in preparation)

The code presented in this repository tests how palaeolithic peoples could have hunted dwarf hippopotamus (<em>Phanourios minor</em> 🦛) and dwarf elephants (<em>Palaeoloxodon cypriotes</em> 🐘) to extinction in the Late Pleistocene.

## <a href="https://github.com/cjabradshaw/CyprusHippoElephant/tree/main/scripts">Scripts</a>
R code by Corey Bradshaw (<a href="http://github.com/cjabradshaw">@cjabradshaw</a>), Frédérik Saltré (<a href="http://github.com/fredsaltre">@fredsaltre</a>), and Salvador Herrando-Pérez

### Cohort-based models
- <code>dwarf hippo & elephant extinction dates.R</code>: estimates Signor-Lipps-corrected extinction window for both megafauna species
- <code>base hippo & elephant model.R</code>: stochastic, age-structured demographic projection models for both megafauna species
- <code>offtake hippo & elephant model.R</code>: stochastic models simulating incrementing offtake rates for both megafauna species (requires running 'base' models first)
- <code>meat equivalents hippo & elephant model.R</code>: stochastic models simulating how incrementing population sizes of humans translates to loss of individuals of both megafauna species (requires running 'base' models first)
- <code>dwarf hippo & elephant model gsa.R</code>: global sensitivity analysis using Latin hypercube sampling of stochastic, age-structured demographic projection models for both megafauna species

### <a href="https://github.com/cjabradshaw/CyprusHippoElephant/tree/main/scripts/source">Source functions</a>
- <code>matrixOperators.r</code>: functions for manipulating matrices for population projections
- <code>qualityRating.r</code>: applies quality rating to radiocarbon age estimates
- <code>endRating.r</code>: chooses final quality rating for radiocarbon age estimates

## <a href="https://github.com/cjabradshaw/CyprusHippoElephant/tree/main/data">Data</a>
- <em>phanourios.txt</em>: radiocarbon dates for <em>Phanourios minor</em> from <a href="http://doi.org/10.1371/journal.pone.0134429">Zazzo et al.</a> (2020: <em>PLoS One</em> 10:0134429)
- <em>palaeoloxodon.txt</em>: radiocarbon dates for <em>Palaeoloxodon cypriotes</em> from Wigand, PE & Simmons, AH (1999: The dating of Akrotiri <em>Aetokremnos</em>, in <a href="https://link.springer.com/book/10.1007/b109876"><em>Faunal Extinction in an Island Society. Pygmy Hippopotamus Hunters of Cyprus</em></a>. AH Simmons (ed). Kluwer Academic Publishers, New York. pp. 193-215)
- <em>qx-Nicolaou.csv</em>: life-table estimates of age-specific survival for <em>Phanourios minor</em> from <a href="http://doi.org/10.1016/j.quaint.2020.09.016">Nicolaou et al.</a> (2020: <em>Quat Int</em> 568:55-64)
- <em>ssdHuman.csv</em>: stable-stage distribution of paleaolithic humans from <a href="http://doi.org/10.21203/rs.3.rs-3468157/v1">Bradshaw et al.</a> (2023: doi:10.21203/rs.3.rs-3468157/v1)

## R libraries
- <a href="https://github.com/FredSaltre/CRIWM/"><code>Rexinct</code></a> (approach described in this <a href="https://doi.org/10.1016/j.quageo.2023.101489">paper</a>): install via Github – <code>devtools::install_github("FredSaltre/CRIWM/Rextinct")</code>. Note: library name has a typo (i.e., 'Rexinct', not 'Rextinct'), so you need to call <code>library(Rexinct)</code> until it is rectified
- <code>dplyr</code>, <code>plotly</code>, <code>ggpubr</code>, <code>truncnorm</code>, <code>doSNOW</code>, <code>iterators</code>, <code>snow</code>, <code>foreach</code>, <code>lhs</code>, <code>data.table</code>, <code>dismo</code>, <code>gbm</code>
<br>
<p>Χάρη στα θρυλικά λουκάνικα της Λυσού</p>
<br>
<p><a href="https://www.ucy.ac.cy"><img align="bottom-left" src="www/UCypruslogo.png" alt="UCyprus logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="http://www.dainst.org"><img align="bottom-left" src="www/DAIlogo.png" alt="DAI logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="30" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp-2.png" alt="GEL logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="https://www.mncn.csic.es"><img align="bottom-left" src="www/CSIClogo.png" alt="Spanish National Research Council logo" height="50" style="margin-top: 20px"></a> &nbsp; <a href="https://www.cut.ac.cy"><img align="bottom-left" src="www/CUTlogoblack.png" alt="CUT logo" height="50" style="margin-top: 20px"></a><a href="https://www.moa.gov.cy/moa/gsd/gsd.nsf/dmlIndex_en/dmlIndex_en"><img align="bottom-left" src="www/CGSlogo.png" alt="CGS logo" height="45" style="margin-top: 20px"></a></p>
